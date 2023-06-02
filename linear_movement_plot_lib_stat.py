"""Static Plot lib module for linear_movement_vibrations"""


# Copyright (C) 2022  Marc Marschall <discrod:MarschallMarc#6420>
#
# This file may be distributed under the terms of the GNU GPLv3 license.

from matplotlib import pyplot as plt
import numpy as np


def plot_frequencies(
    data,
    outfile,
    measurement_parameters,
    gcmd,
    known_causes
):
    plt.ioff()
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    fig.suptitle(
        f"Vibrations while {measurement_parameters.velocity} mm/s linear movement on {measurement_parameters.axis} axis with {measurement_parameters.accel} mm/s^2",
        wrap=True,
    )
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    ax.set_xlabel("frequency in Hz")
    ax.set_ylabel("response")
    ax.set_xlim(data[0][0], measurement_parameters.f_max)
    for length, name in known_causes:
        ax.axvline(x=measurement_parameters.velocity/length, label=name, lw=1, ls='--')
    ax.plot(data[0], data[1], label="x")
    ax.plot(data[0], data[2], label="y")
    ax.plot(data[0], data[3], label="z")
    ax.grid()
    ax.minorticks_on()
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        fancybox=True,
        shadow=False,
        ncol=4,
    )
    # add second abscissa
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.minorticks_on()
    ax2.tick_params(axis="x", direction="out", pad=-1)
    ax2.set_xticklabels(
        [f"{measurement_parameters.velocity/x:.2f}" for x in ax.get_xticks()]
    )
    ax2.set_xlabel("pattern distance in mm")
    fig.tight_layout(pad=0.9)
    plt.savefig(outfile)
    gcmd.respond_info(f"output written to {outfile}")
    plt.close("all")


def plot_relative_power(data, outfile, measurement_parameters, gcmd):
    data = np.array(data)
    plt.ioff()
    fig, ax = plt.subplots()
    fig.suptitle(
        f"Vibration power for axis {measurement_parameters.axis} with accel {measurement_parameters.accel} mm/s^2",
        wrap=True,
    )
    ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    ax.set_xlabel("velocity in mm/s")
    ax.set_ylabel("relative power")
    # heuristics for visibility
    markersize = 5 if len(data) < 200 else 2
    ax.plot(
        data[:, 0],
        data[:, 1:],
        marker="o",
        label=["x component", "y component", "z component"],
        markersize=markersize,
    )
    ax.plot(
        data[:, 0],
        np.sum(data[:, 1:], axis=1),
        marker="o",
        label="total",
        markersize=markersize,
        lw=1,
    )
    ax.grid()
    ax.minorticks_on()
    ax.legend(loc="best", fancybox=True, shadow=False, ncol=1, title="measurement data")
    fig.tight_layout(pad=0.9)

    fig.savefig(outfile)
    gcmd.respond_info(f"output written to {outfile}")
    plt.close("all")


def plot_peak_frequencies(
    data,
    outfile,
    outfilelog,
    measurement_parameters,
    gcmd,
    known_causes, 
):
    plt.ioff()
    fig, ax = plt.subplots()
    velocities, peak_freqs, peak_ffts = zip(*data)
    velocities = np.concatenate(velocities)
    peak_ffts = np.concatenate(peak_ffts)
    peak_freqs = np.concatenate(peak_freqs)
    # peak_ffts = peak_ffts ** (0.8) try out without rescaling
    min_fft = np.amin(peak_ffts)
    normalized_peak_heights = (peak_ffts - min_fft) / (np.amax(peak_ffts) - min_fft)
    peak_height_to_size = 110 * normalized_peak_heights
    scatter = ax.scatter(
        velocities, peak_freqs, c="black", s=peak_height_to_size, marker="o"
    )

    fig.suptitle(
        f"Vibration peak frequencies for axis {measurement_parameters.axis} with accel {measurement_parameters.accel} mm/s^2",
        wrap=True,
    )
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.18, box.width, box.height * 0.85])
    ax.set_xlabel("velocity in mm/s")
    ax.set_ylabel("peak frequency in Hz")
    ax.set_ylim(0, measurement_parameters.f_max)
    ax.minorticks_on()
    for length, name in known_causes:
        ax.plot(velocities, velocities/length, label=length, lw=1)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        fancybox=True,
        shadow=False,
        ncol=3,
    )
    fig.tight_layout(pad=1)
    plt.savefig(outfile)
    gcmd.respond_info(f"output written to {outfile}")
    ax.set_yscale("log")
    plt.axhline(
        y=measurement_parameters.f_max, color="tab:olive", linestyle="--", label="f_max"
    )
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        fancybox=True,
        shadow=False,
        ncol=3,
    )
    ax.set_autoscaley_on(True)
    ax.minorticks_on()
    plt.autoscale(True)
    fig.suptitle(
        f"Vibration peak frequencies for axis {measurement_parameters.axis} with accel {measurement_parameters.accel} mm/s^2, f_max = {measurement_parameters.f_max} Hz",
        wrap=True,
    )
    plt.savefig(outfilelog)
    gcmd.respond_info(f"output written to {outfilelog}")
    plt.close("all")


def plot_frequency_responses_over_velocity(data, outfile, measurement_parameters, gcmd):
    data = np.array(data)
    plt.ioff()
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    fig.suptitle(
        f"Vibration peak frequencies for axis {measurement_parameters.axis} with accel {measurement_parameters.accel} mm/s^2"
    )

    ax.ticklabel_format(style="sci", axis="z", scilimits=(0, 0))
    # fix visuals with simple heuristic
    line_width = 1 if len(data) < 25 else np.exp(-(len(data)) / 40) + 0.2
    for velocity_sample in data[::-1]:
        x = velocity_sample[1]
        y = velocity_sample[2]
        z = velocity_sample[0]
        ax.plot(x, y, zs=z, zdir="y", c="black", lw=line_width)
    ax.set_xlabel("f in Hz")
    ax.set_zlabel("relative response")
    ax.set_ylabel("velocity")
    fig.tight_layout(pad=0.9)
    plt.savefig(outfile)
    gcmd.respond_info(f"output written to {outfile}")
    plt.close("all")
