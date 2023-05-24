"""Static Plot lib module for linear_movement_vibrations"""


# Copyright (C) 2022  Marc Marschall <discrod:MarschallMarc#6420>
#
# This file may be distributed under the terms of the GNU GPLv3 license.

from matplotlib import pyplot as plt
import numpy as np


def plot_frequencies(
    data,
    outfile,
    velocity,
    axis,
    gcmd,
    d=None,
    step_distance=None,
    rotation_distance=None,
    f_max=120,
):
    plt.ioff()
    fig = plt.figure()
    fig.suptitle(f"Vibrations while {velocity}mm/s linear movement on {axis} axis")
    ax = plt.subplot(111)
    box = ax.get_position()
    # shrink and move up to allow legend beneeth
    ax.set_position([box.x0, box.y0 + box.height * 0.18, box.width, box.height * 0.85])
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    ax.set_xlabel("frequency in Hz")
    ax.set_ylabel("response")
    ax.set_xlim(data[0][0], f_max)
    ax.axvline(x=velocity / 2, label="2gt belt pitch", ls="--", color="tab:brown")
    ax.axvline(x=velocity / 1.21, label="2gt belt teeth width", ls="--", color="black")
    ax.axvline(
        x=velocity / 0.80, label="2gt belt valley width", ls="--", color="tab:cyan"
    )
    ax.axvline(
        x=velocity / 0.40, label="2gt belt flat width", ls="--", color="tab:brown"
    )
    if d is not None:
        ax.axvline(
            velocity / (np.pi * d), label="idler rotation", ls="--", color="tab:gray"
        )
    if step_distance is not None:
        ax.axvline(
            velocity / rotation_distance,
            label="pulley rotation",
            ls="--",
            color="tab:olive",
        )
    if rotation_distance is not None:
        ax.axvline(
            velocity * step_distance / rotation_distance,
            label="motor step",
            ls="--",
            color="tab:pink",
        )
    ax.plot(data[0], data[1], label="x")
    ax.plot(data[0], data[2], label="y")
    ax.plot(data[0], data[3], label="z")
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
    ax2.set_position([box.x0, box.y0 + box.height * 0.18, box.width, box.height * 0.85])
    ax2.tick_params(axis="x", direction="in", pad=-15)
    ax2.set_xticklabels([f"{velocity/x:.2f}" for x in ax.get_xticks()])
    ax2.set_xlabel("pattern distance in mm")
    plt.savefig(outfile)
    gcmd.respond_info(f"output written to {outfile}")
    plt.close("all")


def plot_relative_power(data, outfile, axis, gcmd):
    data = np.array(data)
    plt.ioff()
    plt.title(f"Vibration power for axis {axis}")
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    plt.xlabel("velocity in mm/s")
    plt.ylabel("relative power")
    plt.plot(data[:, 0], data[:, 1], marker="o", label="measurement data")
    plt.savefig(outfile)
    gcmd.respond_info(f"output written to {outfile}")
    plt.close("all")


def plot_peak_frequencies(
    data,
    outfile,
    outfilelog,
    axis,
    gcmd,
    d=None,
    step_distance=None,
    rotation_distance=None,
    f_max=200,
):
    plt.ioff()
    fig, ax = plt.subplots()
    velocities, peak_freqs, peak_ffts = zip(*data)
    velocities = np.concatenate(velocities)
    peak_ffts = np.concatenate(peak_ffts)
    peak_freqs = np.concatenate(peak_freqs)
    peak_ffts = peak_ffts ** (0.8)
    normalized_peak_heights = (peak_ffts - np.amin(peak_ffts)) / (
        np.amax(peak_ffts) - np.amin(peak_ffts)
    )
    peak_height_to_size = 90 * normalized_peak_heights
    scatter = ax.scatter(
        velocities, peak_freqs, c="black", s=peak_height_to_size, marker="o"
    )

    fig.suptitle(f"Vibration peak frequencies for axis {axis}")

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.18, box.width, box.height * 0.85])
    ax.set_xlabel("velocity in mm/s")
    ax.set_ylabel("peak frequency in Hz")
    ax.set_ylim(0, f_max)
    ax.plot(velocities, velocities / 2, label="2gt belt pitch")
    ax.plot(velocities, velocities / 1.21, label="2gt belt teeth width")
    ax.plot(velocities, velocities / 0.8, label="2gt belt valley width")
    ax.plot(velocities, velocities / 0.4, label="2gt belt valley flat width")
    if d is not None:
        ax.plot(velocities, velocities / (np.pi * d), label="idler rotation")
    if step_distance is not None:
        ax.plot(velocities, velocities / rotation_distance, label="pulley rotation")
    if rotation_distance is not None:
        ax.plot(
            velocities,
            velocities * step_distance / rotation_distance,
            label="motor step",
        )
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        fancybox=True,
        shadow=False,
        ncol=3,
    )
    plt.savefig(outfile)
    gcmd.respond_info(f"output written to {outfile}")
    ax.set_yscale("log")
    plt.axhline(y=f_max, color="tab:olive", linestyle="--", label="f_max")
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        fancybox=True,
        shadow=False,
        ncol=3,
    )
    ax.set_autoscaley_on(True)
    plt.autoscale(True)
    fig.suptitle(f"Vibration peak frequencies for axis {axis}, f_max = {f_max}Hz  ")
    plt.savefig(outfilelog)
    gcmd.respond_info(f"output written to {outfilelog}")
    plt.close("all")


def plot_frequency_responses_over_velocity(data, outfile, axis, gcmd):
    data = np.array(data)
    plt.ioff()
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    fig.suptitle(f"Vibration peak frequencies for axis {axis}")

    ax.ticklabel_format(style="sci", axis="z", scilimits=(0, 0))

    for velocity_sample in data:
        x = velocity_sample[1]
        y = velocity_sample[2]
        z = velocity_sample[0]
        ax.plot(x, y, zs=z, zdir="y")
    ax.set_xlabel("f in Hz")
    ax.set_zlabel("relative response")
    ax.set_ylabel("velocity")
    plt.savefig(outfile)
    gcmd.respond_info(f"output written to {outfile}")
    plt.close("all")
