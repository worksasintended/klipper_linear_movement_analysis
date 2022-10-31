# vibration measurements for linear, non accelerated movements
# reads adxl and analyses the output for vibrations at linear movements on a defined axis.
# Acceleration phase is not part of the analysis
#
# Copyright (C) 2022  Marc Marschall <discrod:MarschallMarc#6420>
#
# This file may be distributed under the terms of the GNU GPLv3 license.

import datetime
from operator import add

import matplotlib

matplotlib.use("Agg")
import numpy as np
import os
from matplotlib import pyplot as plt


# calculates the mean square of 3d acceleration data summing up all three components
# @param data: array[[t0, x0, y0, z0],...,[tn, xn, yn, zn]]
# @return pd::float
def calculate_total_power(data):
    pd = 0
    norm = len(data)
    for t, x, y, z in data:
        pd += (abs(x) + abs(y) + abs(z)) * (abs(x) + abs(y) + abs(z)) / norm
    return pd


# calculates the frequency spectrum via fft in a given dataset
# @param data: array[[t0, x0, y0, z0],...,[tn, xn, yn, zn]]
# @param f_max::float : maximum frequency considered
def calculate_frequencies(data, f_max):
    frequency_response = []
    start_pos = 0
    end_pos = 0
    dt = (data[len(data) - 1][0] - data[0][0]) / (len(data) - 1)
    for axis in range(1, 4):
        ord_fourier = np.abs(np.fft.rfft(data[:, axis]))
        absc_fourier = np.fft.rfftfreq(data[:, axis].size, dt)
        if start_pos == 0:
            start_pos = np.argmax(absc_fourier > 10)
            end_pos = np.argmax(absc_fourier > f_max)
        frequency_response.append([absc_fourier[start_pos:end_pos], ord_fourier[start_pos:end_pos]])
    return frequency_response


class LinearMovementVibrationsTest:
    def __init__(self, config):
        self.printer = config.get_printer()
        self.printer.register_event_handler("klippy:connect", self.connect)
        self.gcode = self.printer.lookup_object('gcode')
        description = "Measure vibrations at linear movements. Usage:"
        self.gcode.register_command("MEASURE_LINEAR_VIBRATIONS",
                                    self.cmd_MEASURE_LINEAR_VIBRATIONS,
                                    description)
        description = "Measure vibrations at linear movements in a velocity range. Usage:"
        self.gcode.register_command("MEASURE_LINEAR_VIBRATIONS_RANGE",
                                    self.cmd_MEASURE_LINEAR_VIBRATIONS_RANGE,
                                    description)
        # get accel chips, source: resonance_tester.py, should be refactored into helper function on merge
        if not config.get('accel_chip_x', None):
            self.accel_chip_names = [('xy', config.get('accel_chip').strip())]
        else:
            self.accel_chip_names = [
                ('x', config.get('accel_chip_x').strip()),
                ('y', config.get('accel_chip_y').strip())]
            if self.accel_chip_names[0][1] == self.accel_chip_names[1][1]:
                self.accel_chip_names = [('xy', self.accel_chip_names[0][1])]

        self.out_directory = config.get('output_directory')
        self.limits = self._get_limits(config)

    def cmd_MEASURE_LINEAR_VIBRATIONS_RANGE(self, gcmd):
        axis = self._get_axis(gcmd)
        motion_report = self.printer.lookup_object('motion_report')
        v_min, v_max, v_step = self._get_velocity_range(gcmd)
        powers = []
        peak_frequencies = []
        start_pos, end_pos = self._get_move_positions(axis, self.limits)
        for velocity in range(v_min, v_max, v_step):
            measurement_data = self._measure_linear_movement_vibrations(velocity, start_pos, end_pos, motion_report)
            frequency_response = np.array(calculate_frequencies(measurement_data, 1000))
            summed_max_index = self._find_max_total_acceleration(frequency_response)
            peak_frequency = frequency_response[0][0][summed_max_index]
            peak_frequencies.append([velocity, peak_frequency])
            power = calculate_total_power(measurement_data)
            powers.append([velocity, power])
            start_pos_last = start_pos
            start_pos = end_pos
            end_pos = start_pos_last
        if not os.path.exists(self.out_directory):
            os.makedirs(self.out_directory)
        outfile = self._get_outfile_name(self.out_directory, "relative_power")
        self._plot_relative_power(powers, outfile, axis, gcmd)
        outfile = self._get_outfile_name(self.out_directory, "peak_frequencies")
        self._plot_peak_frequencies(peak_frequencies, outfile, axis, gcmd)

    def cmd_MEASURE_LINEAR_VIBRATIONS(self, gcmd):
        axis = self._get_axis(gcmd)
        velocity = self._get_velocity(gcmd)
        motion_report = self.printer.lookup_object('motion_report')
        start_pos, end_pos = self._get_move_positions(axis, self.limits)
        measurement_data = self._measure_linear_movement_vibrations(velocity, start_pos, end_pos, motion_report)
        frequency_response = calculate_frequencies(measurement_data, 200)
        if not os.path.exists(self.out_directory):
            os.makedirs(self.out_directory)
        outfile = self._get_outfile_name(self.out_directory, ("linear_movement_responce_" + str(velocity) + "mmps_"))
        self._plot_frequencies(frequency_response, outfile, velocity, axis, gcmd)

    def _measure_linear_movement_vibrations(self, velocity, start_pos, end_pos, motion_report):
        accel = self.toolhead.max_accel
        self.gcode.run_script_from_command("SET_VELOCITY_LIMIT ACCEL={} ACCEL_TO_DECEL={}".format(accel, accel))
        x_pos, y_pos, z_pos, e_pos = self.toolhead.get_position()
        self.toolhead.move([start_pos[0], start_pos[1], z_pos, e_pos], velocity)
        self.toolhead.wait_moves()
        measurement_handler = [(adxl_axis_attached, accel_chip.start_internal_client())
                               for adxl_axis_attached, accel_chip in self.accel_chips]
        self.toolhead.move([end_pos[0], end_pos[1], z_pos, e_pos], velocity)
        self.toolhead.wait_moves()
        measurement_data = []
        # stop measurement
        for adxl_axis_attached, accel_chip_client in measurement_handler:
            accel_chip_client.finish_measurements()
            if not accel_chip_client.has_valid_samples():
                raise self.gcode.error("No data received from accelerometer")
            else:
                measurement_data = np.asarray(accel_chip_client.get_samples())

        measurement_data_stripped = self._strip_to_linear_velocity_share(velocity, measurement_data, motion_report)
        return measurement_data_stripped

    @staticmethod
    def _strip_to_linear_velocity_share(velocity, data, motion_report):
        # find time stamp of linear movement start
        for i in range(len(data)):
            if motion_report.trapqs['toolhead'].get_trapq_position(data[i, 0])[1] == velocity:
                data = data[i:]
                break
        for i in range(len(data)):
            if motion_report.trapqs['toolhead'].get_trapq_position(data[i, 0])[1] < velocity:
                data = data[0:(i - 1)]
                break
        return data

    @staticmethod
    def _find_max_total_acceleration(frequency_response):
        mapped_frequency_response = map(add, map(add, frequency_response[0][1], frequency_response[1][1]),
                                        frequency_response[2][1])
        summed_max_index = np.argmax(mapped_frequency_response)
        return summed_max_index

    @staticmethod
    def _plot_frequencies(data, outfile, velocity, axis, gcmd):
        plt.ioff()
        plt.title("Vibrations while {}mm/s linear movement on {} axis".format(velocity, axis))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.xlabel("frequency in Hz")
        plt.ylabel("response")
        plt.plot(data[0][0], data[0][1], label="x")
        plt.plot(data[1][0], data[1][1], label="y")
        plt.plot(data[2][0], data[2][1], label="z")
        plt.legend()
        plt.savefig(outfile)
        gcmd.respond_info("output written to {}".format(outfile))
        plt.close('all')

    @staticmethod
    def _plot_relative_power(data, outfile, axis, gcmd):
        data = np.array(data)
        plt.ioff()
        plt.title("Vibration power for axis {}".format(axis))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.xlabel("velocity in mm/s")
        plt.ylabel("relative power")
        plt.plot(data[:, 0], data[:, 1])
        plt.savefig(outfile)
        gcmd.respond_info("output written to {}".format(outfile))
        plt.close('all')

    @staticmethod
    def _plot_peak_frequencies(data, outfile, axis, gcmd):
        data = np.array(data)
        plt.ioff()
        plt.title("Vibration peak frequenices for axis {}".format(axis))
        plt.xlabel("velocity in mm/s")
        plt.ylabel("peak frequency in Hz")
        plt.plot(data[:, 0], data[:, 1], label="measurement data")
        plt.plot(data[:, 0], data[:, 0] / 2, label="belt_teeth_frequency")
        plt.legend()
        plt.savefig(outfile)
        gcmd.respond_info("output written to {}".format(outfile))
        plt.close('all')

    @staticmethod
    def _get_limits(config):
        x_min = int(config.get('x_min'))
        x_max = int(config.get('x_max'))
        y_min = int(config.get('y_min'))
        y_max = int(config.get('y_max'))
        return x_min, x_max, y_min, y_max

    @staticmethod
    def _get_move_positions(axis, limits):
        p1_x = p1_y = p2_x = p2_y = 0
        if axis.lower() == "x":
            p1_x = limits[0]
            p1_y = limits[3] / 2
            p2_x = limits[1]
            p2_y = p1_y
        if axis.lower() == "y":
            p1_x = limits[2] / 2
            p1_y = limits[1]
            p2_x = p1_x
            p2_y = limits[3]
        if axis.lower() == "a":
            p1_x = limits[0]
            p1_y = limits[2]
            p2_x = limits[1]
            p2_y = limits[3]
            p2_x, p2_y = (p1_x, p1_y, p2_x, p2_y)
        if axis.lower() == "b":
            p1_x = limits[1]
            p1_y = limits[2]
            p2_x = limits[0]
            p2_y = limits[3]
            p2_x, p2_y = (p1_x, p1_y, p2_x, p2_y)
        return [p1_x, p1_y], [p2_x, p2_y]

    @staticmethod
    def _get_velocity_range(gcmd):
        vmin = gcmd.get_int("VMIN", None)
        vmin = (vmin, 50)[vmin is None]
        vmax = gcmd.get_int("VMAX", None)
        vmax = (vmax, 300)[vmax is None]
        vstep = gcmd.get_int("STEP", None)
        vstep = (vstep, 10)[vstep is None]
        return vmin, vmax, vstep

    def _get_velocity(self, gcmd):
        velocity = gcmd.get_int("VELOCITY", None)
        velocity = (velocity, 150)[velocity is None]
        if self.toolhead.max_velocity < velocity:
            raise gcmd.error("Requested velocity '{}' succeeds printer limits".format(velocity))
        # TODO check if velocity is reachable
        return velocity

    @staticmethod
    def _get_axis(gcmd):
        axis = gcmd.get("AXIS", None)
        axis = (axis, "x")[axis is None]
        if axis.lower() not in 'xyab':
            raise gcmd.error("Unsupported axis'{}'".format(axis))
        return axis

    @staticmethod
    def _get_outfile_name(directory, filename):
        return directory + filename + datetime.datetime.today().isoformat() + ".png"

    def connect(self):
        self.toolhead = self.printer.lookup_object('toolhead')
        # identical to ResonanceTester.connect, should be moved to helper function on merge
        self.accel_chips = [
            (chip_axis, self.printer.lookup_object(chip_name))
            for chip_axis, chip_name in self.accel_chip_names]

# checks if movement is diagonal and corrects its destination to be a diagonal movement if not
# by using the point towards the center which is the closest to the original fulfilling diagonality
# @param p1_x x-coordinate of starting point
# @param p2_x x-coordinate of the target point
def verify_and_correct_diagonal_move(p1_x, p1_y, p2_x, p2_y):
    if abs(p1_x - p2_x) > abs(p1_y - p2_y):
        p2_x = p1_x + p2_y - p1_y
    elif abs(p1_x - p2_x) < abs(p1_y - p2_y):
        p2_y = p1_y + p2_x - p1_x
    return p2_x, p2_y


def load_config(config):
    return LinearMovementVibrationsTest(config)
