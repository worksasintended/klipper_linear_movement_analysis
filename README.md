# Klipper Linear Movement Vibrations Analysis
This is a klipper extension allowing to measure vibrations on linear movements on a different axis. Unlike with previous solutions, the acceleration phases in the moves are not part of the analyzed data. 

## Usage
This extension adds two new GCODE commands:

### `MEASURE_LINEAR_VIBRATIONS`

**In most usecases it is sufficient to only use `MEASURE_LINEAR_VIBRATIONS [VELOCITY=<velocity>] [AXIS=<x|y|a|b>]`**

This will measure the vibrations frequency spectrum and create a file in the directory defined in the config as follows:

![linear_movement_responce_500mmps_2022-11-06T17_09_49 535678](https://user-images.githubusercontent.com/20718963/200202115-f2bc1d8b-4c0a-4628-9fce-29a9f4677a4b.png)

Full set of options:

`MEASURE_LINEAR_VIBRATIONS [VELOCITY=<velocity>] [AXIS=<x|y|a|b>] [FMAX=<maximum frequency considered default 120>] [FMIN=<minimum frequency considered default 2xVELOCITY>]  [D_IDLER=<diameter of idler>] [XMIN=<VALUE>] [XMAX=<VALUE>] [YMIN=<VALUE>] [YMAX=<VALUE>] [STARTX=<VALUE>] [STARTY=<VALUE>] [ENDX=<VALUE>] [ENDY=<VALUE>] [EXPORT_FFTDATA=<1|0 (enabled|disabled) default is 0>]` 


The settings `XMIN, XMAX, YMINN, YMAX` overwrite the measurement rectangle defined in the config. 

The settings `STARTX, STARTY, ENDX, ENDY`allow to define a movement between any two points. Be aware that this causes the movement to possibly happen on a path different from the defined axis. In this case, the axis is only used to find the adxl corresponding to the axis in the case of multiple accelerometers. Also, make sure that the defined points have a distance big enough to reach the target velocity. 

`D_IDLER` defines the diameter of an idler in mm. If set, the peak frequency graph will show a frequency corresponding to a full rotation of the set idler. This is a way to identify decentered idlers. 
Please be aware, that those frequencies are usually pretty low. To get useful results measure at the fastest speed possible and set `FMIN` to a low value as well as you need to choose the longest travel distance possible.

`FMIN, FMAX` define the frequency range considered. All data outside this range are ignored. Be aware, due to the nature of an fft it does not make sense to use extremely low frequencies that correspond to less than half the measuring time. 

`EXPORT_FFTDATA=1` enables fft data file exportation. In this case, a .npz file in output_directory/raw_data is returned, which can be locally imported after downloading by using e.g. data = numpy.load('fname.npz', allow_pickle=True)['data'].  data array has the form [[velocity, frequencies, fft_data] for velocity in velocities].

### `MEASURE_LINEAR_VIBRATIONS_RANGE`

**A minimal, and in most cases sufficient approach is to use `MEASURE_LINEAR_VIBRATIONS_RANGE [AXIS=<x|y|a|b>] [VMIN=<minimal velocity>] [VMAX=<maximal velocity>] [STEP=<steps size of velocity changes>] `**

This goes through a range of velocities, plots the frequency responses, and calculates the power of the vibrations as well as the frequencies of the main peak for each tested velocity, creating plots as follows:

![frequency_responses_v-range2022-11-06T17_08_02 433594](https://user-images.githubusercontent.com/20718963/200202266-a883232b-1224-411b-a94e-f77ac19949a1.png)
![peak_frequencies_logscale2022-11-06T17_07_59 372651](https://user-images.githubusercontent.com/20718963/200202268-af71abc2-f7da-4b48-abc4-52446ad53799.png)
![peak_frequencies2022-11-06T17_07_59 372548](https://user-images.githubusercontent.com/20718963/200202269-74b2b992-c81d-4a02-8ba2-dcfa5e0c7d72.png)
![relative_power2022-11-06T17_07_58 922808](https://user-images.githubusercontent.com/20718963/200202270-86e9d408-2246-4992-bf54-3dbf3c8bc380.png)


Full set of options:

`MEASURE_LINEAR_VIBRATIONS_RANGE [AXIS=<x|y|a|b>] [VMIN=<minimal velocity>] [VMAX=<maximal velocity>] [STEP=<steps size of veloctity changes>] [D_IDLER=<diameter of pulley or idler>] [FMIN=<minimum frequency considered default 5>] [FMAX=<maximum frequency considered default two times VMAX>] [XMIN=<VALUE>] [XMAX=<VALUE>] [YMIN=<VALUE>] [YMAX=<VALUE>] [STARTX=<VALUE>] [STARTY=<VALUE>] [ENDX=<VALUE>] [ENDY=<VALUE>] [EXPORT_FFTDATA=<1|0 (enabled|disabled) default is 0>]` 


Please read above about the different options, as most of them are identical to `MEASURE_LINEAR_VIBRATIONS`.




## Installation

Log into your raspberry pi via ssh. Clone the git repo via
```
git clone https://github.com/worksasintended/klipper_linear_movement_analysis.git
```
Copy, paste, and excecute the following command:
``` 
bash klipper_linear_movement_analysis/install.sh
```
This will install the extension. The installation of matplotlib and scipy requires some time and creates significant load on the rpi, I suggest NOT doing this, while the rpi is used for printing. 

To update the package via Moonraker or a web frontend, add this to your `moonraker.conf`:
```
[update_manager LinearMovevementAnalysis]
type: git_repo
path: ~/klipper_linear_movement_analysis
primary_branch: development
origin: https://github.com/worksasintended/klipper_linear_movement_analysis.git
install_script: install.sh
env: ~/klippy-env/bin/python
requirements: requirements.txt
managed_services: klipper

```
The moonraker config can either be accessed via ssh or one of the common klipper webfronts, like Fluidd or Mainsail. 


Add the section `[linear_movement_vibrations]` to your `printer.cfg` to activate the plugin
Example configuration:
```
[linear_movement_vibrations]
accel_chip: adxl345
x_min: 5
x_max: 175
y_min: 5
y_max: 175
output_directory: /home/pi/printer_data/config/linear_vibrations/
```
Make sure the defined output directory is writable. The one in the example configuration shown above will create a folder `linear_vibrations` that can be accessed via the file browser in the web frontend, assuming you are using RatOS. Similar to the `input_shaper folder`, you can find it in the machine tab in the `config` root.  If you are not using RatOS and are unsure which directory to use `/tmp/` is a save bet. Be aware, that the pngs will not be automatically removed.
`x_min, x_max,y_min, y_max`define a triangle in which the measurements will be performed. 

If you are using multiple accelerometers, you can also define them as such (untested feature)
```
accel_chip_x: adxl345 rpi
accel_chip_y: adxl345
```
**After these steps, reboot your raspberry, as otherwise the newly installed python libraries are not loaded into the env.**  

## Coding style
- For docstrings, we follow NumPy's "Docstring Standard" and Python's "Docstring Conventions".
- For formatting, we use [`black`](https://black.readthedocs.io/). To implement `black` in your IDE, [find your editor](https://black.readthedocs.io/en/stable/integrations/editors.html)
