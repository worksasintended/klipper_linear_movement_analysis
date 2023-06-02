""" Vibration measurements for linear, non accelerated movements
    reads adxl and analyses the output for vibrations at linear movements 
    on a defined axis. Acceleration phase is not part of the analysis
    """

# Copyright (C) 2022  Marc Marschall <discrod:MarschallMarc#6420>
#
# This file may be distributed under the terms of the GNU GPLv3 license.


import datetime
import os
import enum
import numpy as np
from scipy import signal
from . import linear_movement_plot_lib_stat as plotlib
from dataclasses import dataclass


def calculate_total_power(data):
    """Calculate the mean square of 3d acceleration data for each vibration component.
        The acceleration data is corrected by subtracting the mean (assuming mean = earth accel)

    measurement_parameters
    ----------
    data :  array[[t0, x0, y0, z0],...,[tn, xn, yn, zn]]

    Returns
    -------
    pd : positive float
    """

    norm = len(data)
    mean = np.mean(data, axis=0)
    pd = ((data - mean) ** 2).sum(axis=0) / norm
    return pd


def calculate_frequencies(data, f_max, f_min):
    """Calculates the frequency spectrum via fft in a given dataset.

    measurement_parameters
    ----------
    data :  array[[t0, x0, y0, z0],...,[tn, xn, yn, zn]]
    f_max : float
    f_min : float

    Returns
    -------
    frequency_response : array_like with shape [:, 4]
        Frequency on axis 0, response the axis 1, 2 and 3

    """
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


def verify_and_correct_diagonal_move(p1_x, p1_y, p2_x, p2_y):
    """Correct the destination in case of non-diagonal movement
    by using the point towards the center which is the closest
    to the original fulfilling diagonality.

    measurement_parameters
    ----------
    p1_x p1_y : x and y coordinate start point
    p2_x p2_y : x and y coordinate end point


    Returns
    -------
    p2_x p2_y : adjusted x and y coordinate

    """
    if abs(p1_x - p2_x) > abs(p1_y - p2_y):
        p2_x = p1_x + p2_y - p1_y
    elif abs(p1_x - p2_x) < abs(p1_y - p2_y):
        p2_y = p1_y + p2_x - p1_x
    return p2_x, p2_y


def parse_full_step_distance(config, units_in_radians=None, note_valid=False):
    """source: stepper.py"""

    if units_in_radians is None:
        # Caller doesn't know if units are in radians - infer it
        rd = config.get("rotation_distance", None, note_valid=False)
        gr = config.get("gear_ratio", None, note_valid=False)
        units_in_radians = rd is None and gr is not None
    if units_in_radians:
        rotation_dist = 2.0 * np.pi
        config.get("gear_ratio", note_valid=note_valid)
    else:
        rotation_dist = config.getfloat(
            "rotation_distance", above=0.0, note_valid=note_valid
        )
    # Newer config format with rotation_distance
    full_steps = config.getint(
        "full_steps_per_rotation", 200, minval=1, note_valid=note_valid
    )
    if full_steps % 4:
        raise config.error(
            f"full_steps_per_rotation invalid in section '{config.get_name()}'"
        )
    gearing = parse_gear_ratio(config, note_valid)
    return rotation_dist, full_steps * gearing


def parse_gear_ratio(config, note_valid):
    """source: stepper.py"""

    gear_ratio = config.getlists(
        "gear_ratio", (), seps=(":", ","), count=2, parser=float, note_valid=note_valid
    )
    result = 1.0
    for g1, g2 in gear_ratio:
        result *= g1 / g2
    return result


class GcommandExitType(enum.Enum):
    error = "error"
    success = "success"


class LinearMovementVibrationsTest:
    """This is a klipper extension allowing to measure vibrations on linear movements
    on different axis. Unlike previous solutions, the acceleration phases in the moves
    are not part of the analysed data.

    This extension adds two new GCODE commands:

    - `MEASURE_LINEAR_VIBRATIONS [VELOCITY=<velocity>] [AXIS=<x|y|a|b>]
        [FMAX=<maximum frequency considered default 120>] [FMIN=<minimum frequency
        considered default 2xVELOCITY>]  [D_IDLER=<diameter of idler>] [XMIN=<VALUE>]
        [XMAX=<VALUE>] [YMIN=<VALUE>] [YMAX=<VALUE>] [STARTX=<VALUE>] [STARTY=<VALUE>]
        [ENDX=<VALUE>] [ENDY=<VALUE>] [EXPORT_FFTDATA=<1|0 (enabled|disabled) default is 0>]
        [ACCEL=<set acceleration default max_accel>]`


    - `MEASURE_LINEAR_VIBRATIONS_RANGE [AXIS=<x|y|a|b>] [VMIN=<minimal velocity>]
        [VMAX=<maximal velocity>] [STEP=<steps size of veloctity changes>] [D_IDLER=<diameter
        of pulley or idler>] [FMIN=<minimum frequency considered default 5>] [FMAX=<maximum
        frequency considered default two times VMAX>] [XMIN=<VALUE>] [XMAX=<VALUE>] [YMIN=<VALUE>]
        [YMAX=<VALUE>] [STARTX=<VALUE>] [STARTY=<VALUE>] [ENDX=<VALUE>] [ENDY=<VALUE>]
        [EXPORT_FFTDATA=<1|0 (enabled|disabled) default is 0>]
        [FREQS_PER_V=<number of freqs per velocity> default is 3]
        [ACCEL=<set acceleration default max_accel>]`
    """

    @dataclass()
    class MeasurementParameters:
        axis: str
        v_min: int
        v_max: int
        v_step: int
        velocity: float
        accel: float
        f_max: float
        f_min: float
        start_pos: list
        end_pos: list
        limits: tuple
        freqs_per_v: int

    def __init__(self, config):
        self.toolhead = None
        self.accel_chips = None
        self.max_accel = None
        self.printer = config.get_printer()
        self.printer.register_event_handler("klippy:connect", self.connect)
        self.gcode = self.printer.lookup_object("gcode")
        description = "Measure vibrations at linear movements. Usage:"
        self.gcode.register_command(
            "MEASURE_LINEAR_VIBRATIONS", self.cmd_MEASURE_LINEAR_VIBRATIONS, description
        )
        description = (
            "Measure vibrations at linear movements in a velocity range. Usage:"
        )
        self.gcode.register_command(
            "MEASURE_LINEAR_VIBRATIONS_RANGE",
            self.cmd_MEASURE_LINEAR_VIBRATIONS_RANGE,
            description,
        )
        # get accel chips, source: resonance_tester.py, should be refactored into helper function on merge
        if not config.get("accel_chip_x", None):
            self.accel_chip_names = [("xy", config.get("accel_chip").strip())]
        else:
            self.accel_chip_names = [
                ("x", config.get("accel_chip_x").strip()),
                ("y", config.get("accel_chip_y").strip()),
            ]
            if self.accel_chip_names[0][1] == self.accel_chip_names[1][1]:
                self.accel_chip_names = [("xy", self.accel_chip_names[0][1])]

        self.out_directory = config.get("output_directory")
        if not os.path.exists(self.out_directory):
            os.makedirs(self.out_directory)
        self.limits = self._get_limits_from_config(config)
        self.stepper_configs = self._get_stepper_configs(config)

    def cmd_MEASURE_LINEAR_VIBRATIONS_RANGE(self, gcmd):
        measurement_parameters = self._get_measurement_parameters(gcmd)
        motion_report = self.printer.lookup_object("motion_report")
        velocity_range = range(
            measurement_parameters.v_min,
            measurement_parameters.v_max + 1,
            measurement_parameters.v_step,
        )
        powers = np.zeros((len(velocity_range), 4))

        peak_frequencies = []
        frequency_responses = []

        for vel_idx, velocity in enumerate(velocity_range):
            gcmd.respond_info(f"measuring {velocity} mm/s")
            measurement_parameters.velocity = velocity
            # collect data and add them to the sets
            measurement_data = self._measure_linear_movement_vibrations(
                measurement_parameters, motion_report
            )
            frequency_response = np.array(
                calculate_frequencies(
                    measurement_data,
                    measurement_parameters.f_max,
                    measurement_parameters.f_min,
                )
            )
            mapped_frequency_response = self._map_r3_response_to_single_axis(
                frequency_response
            )
            frequency_responses.append(
                [velocity, frequency_response[0], mapped_frequency_response]
            )

            # peak_properties is a dictionary which contains the peak information (see scipy.signals.find_peaks docu)
            peak_idxs, peak_properties = signal.find_peaks(
                mapped_frequency_response,
                height=(0.0 * np.amax(mapped_frequency_response),),
                distance=1,
            )
            mapped_frequency_response_peaks = mapped_frequency_response[peak_idxs]
            if measurement_parameters.freqs_per_v > len(peak_idxs):
                measurement_parameters.freqs_per_v = -1
            freqs_per_vs = np.argpartition(
                mapped_frequency_response_peaks,
                int(-measurement_parameters.freqs_per_v),
            )[int(-measurement_parameters.freqs_per_v):]
            peak_frequencies.append(
                [
                    np.repeat(velocity, len(freqs_per_vs)),
                    frequency_response[0][peak_idxs][freqs_per_vs],
                    mapped_frequency_response_peaks[freqs_per_vs],
                ]
            )

            power = calculate_total_power(measurement_data)
            powers[vel_idx, 1:4] = power[1:4]
            powers[vel_idx, 0] = velocity
            start_pos_last = measurement_parameters.start_pos
            measurement_parameters.start_pos = measurement_parameters.end_pos
            measurement_parameters.end_pos = start_pos_last

        self._export_fft_data(frequency_responses, gcmd, self.out_directory, "frequency_responses")

        outfile = self._get_outfile_name(self.out_directory, "relative_power")
        plotlib.plot_relative_power(powers, outfile, measurement_parameters, gcmd)
        outfile = self._get_outfile_name(self.out_directory, "peak_frequencies")
        outfilelog = self._get_outfile_name(
            self.out_directory, "peak_frequencies_logscale"
        )
        rotation_dist, step_distance = self._get_step_distance(
            measurement_parameters.axis, self.stepper_configs
        )
        plotlib.plot_peak_frequencies(
            peak_frequencies,
            outfile,
            outfilelog,
            measurement_parameters,
            gcmd,
            d=gcmd.get_float("D_IDLER", None),
            step_distance=step_distance,
            rotation_distance=rotation_dist,
        )
        outfile = self._get_outfile_name(
            self.out_directory, "frequency_responses_v-range"
        )
        plotlib.plot_frequency_responses_over_velocity(
            frequency_responses, outfile, measurement_parameters, gcmd
        )
        self._exit_gcommand()

    def cmd_MEASURE_LINEAR_VIBRATIONS(self, gcmd):
        measurement_parameters = self._get_measurement_parameters(gcmd)
        motion_report = self.printer.lookup_object("motion_report")
        gcmd.respond_info(f"measuring {measurement_parameters.velocity} mm/s")
        measurement_data = self._measure_linear_movement_vibrations(
            measurement_parameters, motion_report
        )

        frequency_response = calculate_frequencies(
            measurement_data, measurement_parameters.f_max, measurement_parameters.f_min
        )

        self._export_fft_data(frequency_response, gcmd, self.out_directory, "frequency_response")

        outfile = self._get_outfile_name(
            self.out_directory,
            (
                    "linear_movement_response_"
                    + str(measurement_parameters.velocity)
                    + "mmps_"
            ),
        )

        rotation_dist, step_distance = self._get_step_distance(
            measurement_parameters.axis, self.stepper_configs
        )
        plotlib.plot_frequencies(
            frequency_response,
            outfile,
            measurement_parameters,
            gcmd,
            d=gcmd.get_float("D_IDLER", None),
            step_distance=step_distance,
            rotation_distance=rotation_dist,
        )
        self._exit_gcommand()

    def _get_init_adxl_handler(self):
        adxl_handler = [
            (adxl_axis_attached, accel_chip.start_internal_client())
            for adxl_axis_attached, accel_chip in self.accel_chips
        ]
        return adxl_handler

    def _measure_linear_movement_vibrations(
            self, measurement_parameters, motion_report
    ):
        self.gcode.run_script_from_command(
            f"SET_VELOCITY_LIMIT ACCEL={measurement_parameters.accel} ACCEL_TO_DECEL={measurement_parameters.accel}"
        )
        x_pos, y_pos, z_pos, e_pos = self.toolhead.get_position()
        self.toolhead.move(
            [
                measurement_parameters.start_pos[0],
                measurement_parameters.start_pos[1],
                z_pos,
                e_pos,
            ],
            measurement_parameters.velocity,
        )
        self.toolhead.wait_moves()
        adxl_handler = self._get_init_adxl_handler()
        self.toolhead.move(
            [
                measurement_parameters.end_pos[0],
                measurement_parameters.end_pos[1],
                z_pos,
                e_pos,
            ],
            measurement_parameters.velocity,
        )
        self.toolhead.wait_moves()
        measurement_data = []
        # stop measurement
        for adxl_axis_attached, accel_chip_client in adxl_handler:
            accel_chip_client.finish_measurements()
            if not accel_chip_client.has_valid_samples():
                self._exit_gcommand(GcommandExitType("error"),"No data received from accelerometer")
            else:
                measurement_data = np.asarray(accel_chip_client.get_samples())
                accel_chip_client.finish_measurements()
        measurement_data_stripped = self._strip_to_linear_velocity_share(
            measurement_parameters.velocity, measurement_data, motion_report, self.gcode
        )
        return measurement_data_stripped

    def connect(self):
        self.toolhead = self.printer.lookup_object("toolhead")
        # identical to ResonanceTester.connect, should be moved to helper function on merge
        self.accel_chips = [
            (chip_axis, self.printer.lookup_object(chip_name))
            for chip_axis, chip_name in self.accel_chip_names
        ]
        self.max_accel = self.toolhead.max_accel

    def _exit_gcommand(self, state=GcommandExitType("success"), message=None):
        self.toolhead.max_accel = self.max_accel
        for adxl_axis_attached, accel_chip in self.accel_chips:
            if accel_chip.is_measuring():
                # no way to reach it without using protected function as of now
                accel_chip._finish_measurements()

        if state.value == "error":
            raise self.gcode.error(message)


    def _export_fft_data(self, frequency_response, gcmd, out_directory, fname):
        if gcmd.get_int("EXPORT_FFTDATA", 0) == 1:
            outfile = self._get_outfile_name("", fname, "")
            self._write_data_outfile(
                out_directory, gcmd, outfile, frequency_response
            )
    def _get_accel(self, gcmd, max_accel):
        # define max_accel from toolhead and check if user settings exceed max accel
        accel = gcmd.get_int("ACCEL", max_accel)
        if accel > max_accel:
            accel = max_accel
            gcmd.respond_info(
                f"Warning: Cannot exceed machine limits. Acceleration set to {max_accel} mm/s^2"
            )
        return accel

    def _get_measurement_parameters(self, gcmd):
        axis = self._get_axis(gcmd)
        v_min, v_max, v_step = self._get_velocity_range(gcmd)
        velocity = self._get_velocity(gcmd)
        accel = self._get_accel(gcmd, self.max_accel)
        f_max = gcmd.get_int("FMAX", 2 * v_max)
        f_min = gcmd.get_int("FMIN", 5)
        limits = self._get_limits_from_gcode(gcmd, self.limits)
        start_pos, end_pos = self._get_move_positions(axis, limits, gcmd)
        freqs_per_v = self._get_freqs_per_v(gcmd)

        return self.MeasurementParameters(
            axis,
            v_min,
            v_max,
            v_step,
            velocity,
            accel,
            f_max,
            f_min,
            start_pos,
            end_pos,
            limits,
            freqs_per_v,
        )
    def _strip_to_linear_velocity_share(self,velocity, data, motion_report, gcmd):
        # find time stamp of linear movement start
        velocity_not_reached = True
        for i in range(len(data)):
            if (
                    motion_report.trapqs["toolhead"].get_trapq_position(data[i, 0])[1]
                    == velocity
            ):
                data = data[i:]
                velocity_not_reached = False
                break
        for i in range(len(data)):
            if (
                    motion_report.trapqs["toolhead"].get_trapq_position(data[i, 0])[1]
                    < velocity
            ):
                data = data[: i - 1]
                break
        if velocity_not_reached or len(data) < 300:
            message = "Target velocity not reached for a sufficient amount of time. Either decrease target velocity, " \
                      "increase acceleration or increase test area "
            self._exit_gcommand(GcommandExitType("error"), message)

        return data

    @staticmethod
    def _get_freqs_per_v(gcmd):
        freqs_per_v = gcmd.get_int("FREQS_PER_V", 3)
        if freqs_per_v == 0:
            freqs_per_v = 1
        return freqs_per_v

    @staticmethod
    def _write_data_outfile(directory, gcmd, fname, data):
        """Write data into out_directory/raw_data/fname by np.savez."""

        if not os.path.exists(directory + "raw_data"):
            os.makedirs(directory + "raw_data")
        outfile = directory + "raw_data/" + fname
        np.savez(outfile, data=np.array(data, dtype=object))
        gcmd.respond_info(f"data output written to {outfile}")

    @staticmethod
    def _get_stepper_configs(config):
        stepper_config = []
        for stepper in ["stepper_x", "stepper_y"]:
            stepper_config.append(config.getsection(stepper))
        return stepper_config

    @staticmethod
    def _get_step_distance(axis, config):
        rotation_dist = step_distance = None
        if axis.lower() in "x":
            rotation_dist, step_distance = parse_full_step_distance(config[0])
        elif axis.lower() in "y":
            rotation_dist, step_distance = parse_full_step_distance(config[1])
        return rotation_dist, step_distance




    @staticmethod
    def _map_r3_response_to_single_axis(frequency_response):
        combined_array = np.array(
            [
                frequency_response[1] ** 2,
                frequency_response[2] ** 2,
                frequency_response[3] ** 2,
            ]
        )
        mapped_frequency_response = np.sqrt(combined_array.sum(axis=0))
        return mapped_frequency_response

    @staticmethod
    def _get_limits_from_config(config):
        x_min = int(config.get("x_min"))
        x_max = int(config.get("x_max"))
        y_min = int(config.get("y_min"))
        y_max = int(config.get("y_max"))
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
            message = f"Requested velocity '{velocity}' succeeds printer limits"
            self._exit_gcommand(GcommandExitType("error"), message)
        return velocity


    def _get_axis(self,gcmd):
        axis = gcmd.get("AXIS", None)
        axis = (axis, "x")[axis is None]
        if axis.lower() not in ["x", "y", "a", "b"]:
            self._exit_gcommand(GcommandExitType("error"), f"Unsupported axis'{axis}'")
        return axis

    @staticmethod
    def _get_outfile_name(directory, fname, extension=".png"):
        return directory + fname + "_" + datetime.datetime.today().isoformat() + extension


def load_config(config):
    return LinearMovementVibrationsTest(config)
