# vibration measurements for linear, non accelerated movements
# reads adxl and analyses the output for vibrations at linear movements on a defined axis.
# Acceleration phase is not part of the analysis
#
# Copyright (C) 2022  Marc Marschall <discrod:MarschallMarc#6420>
#
# This file may be distributed under the terms of the GNU GPLv3 license.

import datetime

import matplotlib

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
def calculate_frequencies(data, f_max, f_min):
    dt = (data[len(data) - 1][0] - data[0][0]) / (len(data) - 1)
    norm = data[:, 0].size
    absc_fourier = np.fft.rfftfreq(norm, dt)
    start_pos = np.argmax(absc_fourier >= f_min)
    end_pos = np.argmax(absc_fourier >= f_max)
    if end_pos == 0:
        end_pos = data[:, 0].size
    frequency_response = [absc_fourier[start_pos:end_pos]]
    for axis in range(1, 4):
        ord_fourier = np.abs(np.fft.rfft(data[:, axis]))
        frequency_response.append(ord_fourier[start_pos:end_pos])
    return frequency_response


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


# source: stepper.py
def parse_full_step_distance(config, units_in_radians=None, note_valid=False):
    if units_in_radians is None:
        # Caller doesn't know if units are in radians - infer it
        rd = config.get('rotation_distance', None, note_valid=False)
        gr = config.get('gear_ratio', None, note_valid=False)
        units_in_radians = rd is None and gr is not None
    if units_in_radians:
        rotation_dist = 2. * np.pi
        config.get('gear_ratio', note_valid=note_valid)
    else:
        rotation_dist = config.getfloat('rotation_distance', above=0.,
                                        note_valid=note_valid)
    # Newer config format with rotation_distance
    full_steps = config.getint('full_steps_per_rotation', 200, minval=1,
                               note_valid=note_valid)
    if full_steps % 4:
        raise config.error("full_steps_per_rotation invalid in section '%s'"
                           % (config.get_name(),))
    gearing = parse_gear_ratio(config, note_valid)
    return rotation_dist, full_steps * gearing


# source: stepper.py
def parse_gear_ratio(config, note_valid):
    gear_ratio = config.getlists('gear_ratio', (), seps=(':', ','), count=2,
                                 parser=float, note_valid=note_valid)
    result = 1.
    for g1, g2 in gear_ratio:
        result *= g1 / g2
    return result


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
        self.limits = self._get_limits_from_config(config)
        self.stepper_configs = self._get_stepper_configs(config)

    def cmd_MEASURE_LINEAR_VIBRATIONS_RANGE(self, gcmd):
        axis = self._get_axis(gcmd)
        motion_report = self.printer.lookup_object('motion_report')
        v_min, v_max, v_step = self._get_velocity_range(gcmd)
        f_max = gcmd.get_int("FMAX", 2*v_max)
        powers = []
        peak_frequencies = []
        frequency_responses = []
        limits = self._get_limits_from_gcode(gcmd, self.limits)
        start_pos, end_pos = self._get_move_positions(axis, limits, gcmd)
        for velocity in range(v_min, v_max + 1, v_step):
            gcmd.respond_info("measuring {} mm/s".format(velocity))
            # collect data and add them to the sets
            measurement_data = self._measure_linear_movement_vibrations(velocity, start_pos, end_pos, motion_report)
            frequency_response = np.array(
                calculate_frequencies(measurement_data, f_max, gcmd.get_int("FMIN", 5)))
            mapped_frequency_response = self._map_r3_response_to_single_axis(frequency_response)
            frequency_responses.append([velocity, frequency_response[0], mapped_frequency_response])
            summed_max_index = np.argmax(mapped_frequency_response)
            peak_frequency = frequency_response[0][summed_max_index]
            peak_frequencies.append([velocity, peak_frequency])
            power = calculate_total_power(measurement_data)
            powers.append([velocity, power])
            # reverse movement
            start_pos_last = start_pos
            start_pos = end_pos
            end_pos = start_pos_last
        if not os.path.exists(self.out_directory):
            os.makedirs(self.out_directory)
        outfile = self._get_outfile_name(self.out_directory, "relative_power")
        self._plot_relative_power(powers, outfile, axis, gcmd)
        outfile = self._get_outfile_name(self.out_directory, "peak_frequencies")
        outfilelog = self._get_outfile_name(self.out_directory, "peak_frequencies_logscale")
        rotation_dist, step_distance = self._get_step_distance(axis, self.stepper_configs)
        self._plot_peak_frequencies(peak_frequencies, outfile, outfilelog, axis, gcmd,
                                    d=gcmd.get_float("D_IDLER", None),
                                    step_distance=step_distance, rotation_distance=rotation_dist, f_max=f_max)
        outfile = self._get_outfile_name(self.out_directory, "frequency_responses_v-range")
        self._plot_frequency_responses_over_velocity(frequency_responses, outfile, axis, gcmd)

    def cmd_MEASURE_LINEAR_VIBRATIONS(self, gcmd):
        axis = self._get_axis(gcmd)
        velocity = self._get_velocity(gcmd)
        motion_report = self.printer.lookup_object('motion_report')
        limits = self._get_limits_from_gcode(gcmd, self.limits)
        start_pos, end_pos = self._get_move_positions(axis, limits, gcmd)
        measurement_data = self._measure_linear_movement_vibrations(velocity, start_pos, end_pos, motion_report)
        f_max = gcmd.get_int("FMAX", 2*velocity)
        frequency_response = calculate_frequencies(measurement_data, f_max,
                                                   gcmd.get_int("FMIN", 5))
        if not os.path.exists(self.out_directory):
            os.makedirs(self.out_directory)
        outfile = self._get_outfile_name(self.out_directory, ("linear_movement_responce_" + str(velocity) + "mmps_"))
        rotation_dist, step_distance = self._get_step_distance(axis, self.stepper_configs)
        self._plot_frequencies(frequency_response, outfile, velocity, axis, gcmd, d=gcmd.get_float("D_IDLER", None),
                               step_distance=step_distance, rotation_distance=rotation_dist, f_max=f_max)

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

        measurement_data_stripped = self._strip_to_linear_velocity_share(velocity, measurement_data, motion_report,
                                                                         self.gcode)
        return measurement_data_stripped

    @staticmethod
    def _get_stepper_configs(config):
        stepper_config = []
        for stepper in ['stepper_x', 'stepper_y']:
            stepper_config.append(config.getsection(stepper))
        return stepper_config

    @staticmethod
    def _get_step_distance(axis, config):
        rotation_dist = step_distance = None
        if axis.lower() in 'x':
            rotation_dist, step_distance = parse_full_step_distance(config[0])
        elif axis.lower() in 'y':
            rotation_dist, step_distance = parse_full_step_distance(config[1])
        return rotation_dist, step_distance

    @staticmethod
    def _strip_to_linear_velocity_share(velocity, data, motion_report, gcmd):
        # find time stamp of linear movement start
        velocity_not_reached = True
        for i in range(len(data)):
            if motion_report.trapqs['toolhead'].get_trapq_position(data[i, 0])[1] == velocity:
                data = data[i:]
                velocity_not_reached = False
                break
        for i in range(len(data)):
            if motion_report.trapqs['toolhead'].get_trapq_position(data[i, 0])[1] < velocity:
                data = data[0:(i - 1)]
                break
        if velocity_not_reached or len(data) < 300:
            raise gcmd.error("Target velocity not reached for a sufficient amount of time. Either decrease target "
                             "velocity, increase acceleration or increase test area ")
        return data

    @staticmethod
    def _map_r3_response_to_single_axis(frequency_response):
        combined_array = np.array([frequency_response[1], frequency_response[2], frequency_response[3]])
        mapped_frequency_response = combined_array.sum(axis=0)
        return mapped_frequency_response

    @staticmethod
    def _get_limits_from_config(config):
        x_min = int(config.get('x_min'))
        x_max = int(config.get('x_max'))
        y_min = int(config.get('y_min'))
        y_max = int(config.get('y_max'))
        return x_min, x_max, y_min, y_max

    @staticmethod
    def _get_limits_from_gcode(gcmd, limits):
        x_min = gcmd.get_int("XMIN", limits[0])
        x_max = gcmd.get_int("XMAX", limits[1])
        y_min = gcmd.get_int("YMIN", limits[2])
        y_max = gcmd.get_int("YMAX", limits[3])
        return x_min, x_max, y_min, y_max

    @staticmethod
    def _get_move_positions(axis, limits, gcmd):
        p1_x = p1_y = p2_x = p2_y = 0
        if axis.lower() == "x":
            p1_x = limits[0]
            p1_y = limits[3] / 2
            p2_x = limits[1]
            p2_y = p1_y
        elif axis.lower() == "y":
            p1_x = limits[1] / 2
            p1_y = limits[2]
            p2_x = p1_x
            p2_y = limits[3]
        elif axis.lower() == "a":
            p1_x = limits[0]
            p1_y = limits[2]
            p2_x = limits[1]
            p2_y = limits[3]
            p2_x, p2_y = verify_and_correct_diagonal_move(p1_x, p1_y, p2_x, p2_y)
        elif axis.lower() == "b":
            p1_x = limits[1]
            p1_y = limits[2]
            p2_x = limits[0]
            p2_y = limits[3]
            p2_x, p2_y = verify_and_correct_diagonal_move(p1_x, p1_y, p2_x, p2_y)

        p1_x = gcmd.get_int("STARTX", p1_x)
        p1_y = gcmd.get_int("STARTY", p1_y)
        p2_x = gcmd.get_int("ENDX", p2_x)
        p2_y = gcmd.get_int("ENDY", p2_y)
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

    @staticmethod
    def _plot_frequencies(data, outfile, velocity, axis, gcmd, d=None, step_distance=None, rotation_distance=None,
                          f_max=120):
        plt.ioff()
        fig = plt.figure()
        fig.suptitle("Vibrations while {}mm/s linear movement on {} axis".format(velocity, axis))
        ax = plt.subplot(111)
        box = ax.get_position()
        # shrink and move up to allow legend beneeth
        ax.set_position([box.x0, box.y0 + box.height * 0.18, box.width, box.height * 0.85])
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax.set_xlabel("frequency in Hz")
        ax.set_ylabel("response")
        ax.set_xlim(data[0][0],f_max)
        ax.axvline(x=velocity / 2, label="2gt belt pitch", ls='--', color='tab:brown')
        ax.axvline(x=velocity / 1.21, label="2gt belt teeth width", ls='--', color='black')
        ax.axvline(x=velocity / 0.80, label="2gt belt valley width", ls='--', color='tab:cyan')
        ax.axvline(x=velocity / 0.40, label="2gt belt flat width", ls='--', color='tab:brown')
        if d is not None:
            ax.axvline(velocity / (np.pi * d), label="idler rotation", ls='--', color='tab:gray')
        if step_distance is not None:
            ax.axvline(velocity / rotation_distance, label="pulley rotation", ls='--', color='tab:olive')
        if rotation_distance is not None:
            ax.axvline(velocity * step_distance / rotation_distance, label="motor step", ls='--', color='tab:pink')
        ax.plot(data[0], data[1], label="x")
        ax.plot(data[0], data[2], label="y")
        ax.plot(data[0], data[3], label="z")
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), fancybox=True, shadow=False, ncol=4)
        #add second abscissa
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_position([box.x0, box.y0 + box.height * 0.18, box.width, box.height * 0.85])
        ax2.tick_params(axis="x",direction="in", pad=-15)
        ax2.set_xticklabels(['{0:.2f}'.format(velocity/x) for x in ax.get_xticks()])
        ax2.set_xlabel("pattern distance in mm")
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
        plt.plot(data[:, 0], data[:, 1], marker='o', label="measurement data")
        plt.savefig(outfile)
        gcmd.respond_info("output written to {}".format(outfile))
        plt.close('all')

    @staticmethod
    def _plot_peak_frequencies(data, outfile, outfilelog, axis, gcmd, d=None, step_distance=None,
                               rotation_distance=None, f_max=200):
        data = np.array(data)
        plt.ioff()
        fig = plt.figure()
        fig.suptitle("Vibration peak frequencies for axis {}".format(axis))
        ax = plt.subplot(111)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.18, box.width, box.height * 0.85])
        ax.set_xlabel("velocity in mm/s")
        ax.set_ylabel("peak frequency in Hz")
        ax.set_ylim(0, f_max)
        ax.plot(data[:, 0], data[:, 1], linestyle='--', marker='o', label="measurement data")
        ax.plot(data[:, 0], data[:, 0] / 2, label="2gt belt pitch")
        ax.plot(data[:, 0], data[:, 0] / 1.21, label="2gt belt teeth width")
        ax.plot(data[:, 0], data[:, 0] / .8, label="2gt belt valley width")
        ax.plot(data[:, 0], data[:, 0] / .4, label="2gt belt valley flat width")
        if d is not None:
            ax.plot(data[:, 0], data[:, 0] / (np.pi * d), label="idler rotation")
        if step_distance is not None:
            ax.plot(data[:, 0], data[:, 0] / rotation_distance, label="pulley rotation")
        if rotation_distance is not None:
            ax.plot(data[:, 0], data[:, 0] * step_distance / rotation_distance, label="motor step")
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), fancybox=True, shadow=False, ncol=3)
        plt.savefig(outfile)
        gcmd.respond_info("output written to {}".format(outfile))
        ax.set_yscale('log')
        plt.axhline(y=f_max, color='tab:olive', linestyle='--', label="f_max")
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), fancybox=True, shadow=False, ncol=3)
        ax.set_autoscaley_on(True)
        plt.autoscale(True)
        fig.suptitle("Vibration peak frequencies for axis {}, f_max = {}Hz  ".format(axis, f_max))
        plt.savefig(outfilelog)
        gcmd.respond_info("output written to {}".format(outfilelog))
        plt.close('all')

    @staticmethod
    def _plot_frequency_responses_over_velocity(data, outfile, axis, gcmd):
        data = np.array(data)
        plt.ioff()
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        fig.suptitle("Vibration peak frequencies for axis {}".format(axis))
        ax.ticklabel_format(style='sci', axis='z', scilimits=(0, 0))

        for velocity_sample in data:
            x = velocity_sample[1]
            y = velocity_sample[2]
            z = velocity_sample[0]
            ax.plot(x, y, zs=z, zdir='y')
        ax.set_xlabel("f in Hz")
        ax.set_zlabel("relative response")
        ax.set_ylabel("velocity")
        plt.savefig(outfile)
        gcmd.respond_info("output written to {}".format(outfile))
        plt.close('all')


def load_config(config):
    return LinearMovementVibrationsTest(config)
