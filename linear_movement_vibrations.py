# vibration measurements for linear, non accelerated movements 
#
# Copyright (C) 2022  Marc Marschall <MarschallMarc#6420>
#
# This file may be distributed under the terms of the GNU GPLv3 license.

import numpy as np
import matplotlib, os
matplotlib.use("Agg")

from matplotlib import pyplot as plt
from datetime import datetime   
from operator import add

# TODO:
# error handling
# range handling
# axis handling
# refactoring
# limit handling
# check if velocity is reachable long enough
# define long enough (2cm?)
# power over velocity check

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

    def cmd_MEASURE_LINEAR_VIBRATIONS_RANGE(self,gcmd):
        axis = self._get_axis(gcmd)
        motion_report = self.printer.lookup_object('motion_report')
        v_min, v_max, v_step = self._get_velocity_range(gcmd)
        power_densities = []
        peak_frequencies = []
        for velocity in range(v_min, v_max, v_step):
            measurement_data = self._measure_linear_movement_vibrations(velocity,axis, motion_report)
            # measurement step time
            dt = (measurement_data[len(measurement_data)-1][0]-measurement_data[0][0])/(len(measurement_data)-1)
            frequency_response = np.array(self._calculate_frequencies(measurement_data, dt,1000))
            mapped_frequency_response = map(add, map(add, frequency_response[0][1],frequency_response[1][1]),frequency_response[2][1])
            summed_max_index = np.argmax(mapped_frequency_response)
            peak_frequency = frequency_response[0][0][summed_max_index]
            peak_frequencies.append([velocity,peak_frequency])
            power_density = self._calculate_total_power_density(measurement_data)
            power_densities.append([velocity,power_density])
        out_directory = "/home/pi/klipper_config/linear_vibrations/"
        if not os.path.exists(out_directory):
            os.makedirs(out_directory)
        outfile = out_directory + "power_density" + datetime.today().isoformat() + ".png"
        self._plot_power_density(power_densities, outfile, axis)
        outfile = out_directory + "peak_frequencies" + datetime.today().isoformat() + ".png"
        self._plot_peak_frequencies(peak_frequencies, outfile, axis)

    def _get_velocity_range(self,gcmd):
        vmin = int(gcmd.get("VMIN", None))
        vmin = (vmin, 100)[vmin is None]
        vmax = int(gcmd.get("VMAX", None))
        vmax = (vmax, 200)[vmax is None]
        vstep = int(gcmd.get("STEP", None))
        vstep = (vstep, 10)[vstep is None]
        return vmin, vmax, vstep
    
    def _plot_power_density(self, data, outfile,axis):
        data = np.array(data)
        plt.ioff()
        plt.title("Vibration power density for axis {}".format(axis))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.xlabel("velocity in mm/s")
        plt.ylabel("power density")
        plt.plot(data[:,0],data[:,1])
        plt.savefig(outfile)
        self.gcode.respond_info("output written to {}".format(outfile))
        plt.close('all')
        
    def _plot_peak_frequencies(self, data, outfile,axis):
        data = np.array(data)
        plt.ioff()
        plt.title("Vibration peak frequenices for axis {}".format(axis))
        plt.xlabel("velocity in mm/s")
        plt.ylabel("peak frequency in Hz")
        plt.plot(data[:,0],data[:,1], label="measurement data")
        plt.plot(data[:,0], data[:,0]/2, label="belt_teeth_frequency")
        plt.legend()
        plt.savefig(outfile)
        self.gcode.respond_info("output written to {}".format(outfile))
        plt.close('all')

    def _calculate_total_power_density(self, data):
        pd = 0
        norm = len(data)
        for t, x, y, z in data:
            pd+=abs(x)/norm
            pd+=abs(y)/norm
            pd+=abs(z)/norm
        return pd

    def cmd_PRINT_LAST_MOVE_TIME(self,gcmd):
        self.gcode.respond_info(str(self.toolhead.get_last_move_time()))

    def cmd_MEASURE_LINEAR_VIBRATIONS(self, gcmd):
        axis = self._get_axis(gcmd)
        velocity = self._get_velocity(gcmd)
        motion_report = self.printer.lookup_object('motion_report')
        measurement_data = self._measure_linear_movement_vibrations(velocity,axis, motion_report)
        dt = (measurement_data[len(measurement_data)-1][0]-measurement_data[0][0])/(len(measurement_data)-1)
        #TODO manage axis assignments
        frequency_response = self._calculate_frequencies(measurement_data, dt,200)
        out_directory = "/home/pi/klipper_config/linear_vibrations/"
        if not os.path.exists(out_directory):
            os.makedirs(out_directory)
        outfile = out_directory + "linear_movement_responce_" + str(velocity) +  "mmps_" + datetime.today().isoformat() + ".png"
        self._plot_frequencies(frequency_response, outfile, velocity, axis)
            
    def _get_velocity(self, gcmd):
        velocity = int(gcmd.get("VELOCITY", None))
        velocity = (velocity,150)[velocity is None]
        if self.toolhead.max_velocity < velocity:
            raise gcmd.error("Requested velocity '{}' suceeds printer limits".format(velocity))
        #TODO check if velocity is reachable
        return velocity

    def _get_axis(self, gcmd):
        axis = gcmd.get("AXIS", None)
        axis = (axis, "x")[axis is None]
        if axis.lower() not in 'xy':
            raise gcmd.error("Unsupported axis'{}'".format(axis))
        return axis

    def _plot_frequencies(self, data, outfile, velocity, axis):      
        plt.ioff()
        plt.title("Vibrations while {}mm/s linear movement on {} axis".format(velocity,axis))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.xlabel("frequency in Hz")
        plt.ylabel("response")
        plt.plot(data[0][0], data[0][1], label="x")
        plt.plot(data[1][0], data[1][1], label="y")
        plt.plot(data[2][0], data[2][1], label="z")
        plt.legend()
        plt.savefig(outfile)
        self.gcode.respond_info("output written to {}".format(outfile))
        plt.close('all')

    def _calculate_frequencies(self, data, dt, f_max):
        frequency_response = []
        start_pos = 0
        end_pos = 0
        for axis in range(1,4):
            ord_fourier = np.abs(np.fft.rfft(data[:,axis]))
            absc_fourier =  np.fft.rfftfreq(data[:,axis].size, dt)
            if start_pos == 0:
                    start_pos = np.argmax(absc_fourier > 10)
                    end_pos = np.argmax(absc_fourier > f_max)
            frequency_response.append([absc_fourier[start_pos:end_pos],ord_fourier[start_pos:end_pos]])
        return frequency_response

    def _measure_linear_movement_vibrations(self, velocity, axis, motion_report):
        accel = self.toolhead.max_accel
        self.gcode.run_script_from_command("SET_VELOCITY_LIMIT ACCEL={} ACCEL_TO_DECEL={}".format(accel,accel))
        X, Y, Z, E = self.toolhead.get_position()
        xlim, ylim = self._get_limits() 

        measurement_data = []
                

        # move to starting position
        self.toolhead.move([5, ylim/2, Z, E], velocity)
        self.toolhead.wait_moves()

        # start measuremnt
        measurement_handler = [(adxl_axis_attached, accel_chip.start_internal_client())
                      for adxl_axis_attached, accel_chip in self.accel_chips]


        self.toolhead.move([xlim-5, ylim/2, Z, E], velocity)
        self.toolhead.wait_moves()
        measurement_data = []

        #stop measurement
        for adxl_axis_attached, accel_chip_client in measurement_handler:
            accel_chip_client.finish_measurements()
            if not accel_chip_client.has_valid_samples():
                raise self.gcode.error("No data received from acceleratometer")
            else:
                measurement_data = np.asarray(accel_chip_client.get_samples())
        
        measurement_data_stripped = self._strip_to_linear_velocity_share(velocity, measurement_data, motion_report)
        return measurement_data_stripped

    def _get_limits(self): #TODO
        return 180,180

    #TODO
    def _get_direction(self,axis): 
        return 1,0

    def _get_starting_position(self,axis,limits): #TODO
        return 1,1

    def _strip_to_linear_velocity_share(self, velocity, data, motion_report):
        # find time stamp of linear movement start
        for i in range(len(data)):
            if motion_report.trapqs['toolhead'].get_trapq_position(data[i,0])[1] == velocity:
                data = data[i:]
                break
        for i in range(len(data)):
            if motion_report.trapqs['toolhead'].get_trapq_position(data[i,0])[1] < velocity:
                data = data[0:(i-1)]
                break
        return data


    def connect(self):
        self.toolhead = self.printer.lookup_object('toolhead')
        # identical to ResonanceTester.connect, should be moved to helper function on merge
        self.accel_chips = [
                (chip_axis, self.printer.lookup_object(chip_name))
                for chip_axis, chip_name in self.accel_chip_names]


def load_config(config):
    return LinearMovementVibrationsTest(config)
