# Klipper Linear Movement Vibrations Analysis
This is a klipper extension allowing to measure vibrations on linear movements on different axis. Unlike with previous solutions, the acceleration phases in the moves are not part of the analysed data. 

## Usage
This extension adds two new GCODE commands:

`MEASURE_LINEAR_VIBRATIONS [VELOCITY=<velocity>] [AXIS=<x|y|a|b>] [FMAX=<maximum frequency considered default 120>] [FMIN=<minimum frequency considered default 10>] [XMIN=<VALUE>] [XMAX=<VALUE>] [YMIN=<VALUE>] [YMAX=<VALUE>] [STARTX=<VALUE>] [STARTY=<VALUE>] [ENDX=<VALUE>] [ENDY=<VALUE>]` will measure the vibrations frequency spectrum and create a file in the directory defined in the config as follows:

![linear_movement_responce_150mmps_2022-10-30T17_53_59 439905](https://user-images.githubusercontent.com/20718963/199113335-7f21d635-22e4-4c77-abc3-ec5677382064.png)

The settings `XMIN, XMAX, YMINN, YMAX` overwrite the measurement rectangle defined in the config. 

The settings `STARTX, STARTY, ENDX, ENDY`allow to define a movement between any two points. Be aware that this causes the movement to possibly happen on a path different to the defined axis. In this case the axis is only used to find the adxl corresponding to the axis in case of multiple accelerometers. Also make sure that the defined points have a distance big enough to reach target velocity. 

**In most usecases it is sufficient to only use `MEASURE_LINEAR_VIBRATIONS [VELOCITY=<velocity>] [AXIS=<x|y|a|b>]`**

`MEASURE_LINEAR_VIBRATIONS_RANGE [AXIS=<x|y|a|b>] [VMIN=<minimal velocity>] [VMAX=<maximal velocity>] [STEP=<steps size of veloctity changes>] [DROT=<diameter of pulley or idler>] [FMIN=<minimum frequency considered default 10>] [FMAX=<maximum frequency considered default 120>] [XMIN=<VALUE>] [XMAX=<VALUE>] [YMIN=<VALUE>] [YMAX=<VALUE>] [STARTX=<VALUE>] [STARTY=<VALUE>] [ENDX=<VALUE>] [ENDY=<VALUE>]` goes through a range of velocities, plots the frequency responses and calculates the power of the vibrations as well as the frequencies of the main peak for each tested velocity, creating plots as following:

![frequency_responses_v-range2022-11-01T13_55_39 067495](https://user-images.githubusercontent.com/20718963/199251639-0972baed-a081-4a83-aa8b-13150158ad59.png)
![peak_frequencies2022-11-01T14_01_02 980854](https://user-images.githubusercontent.com/20718963/199255538-db5a9a7b-c424-44b3-b473-598dd4df22c5.png)
![relative_power2022-10-31T21_32_58 727183](https://user-images.githubusercontent.com/20718963/199114782-a5c26bd9-f85c-4b45-90e0-74596a00c371.png)

Please read above about the different settings possible. `DROT` defines the diameter of an idler or pulley. If set, the peak frequency graph will show a frequency corresponding to a full rotation of set idler or pulley. This is a way to identify decentered idlers or pulleys.   

**A minimal, and in most cases sufficient approach is to use `MEASURE_LINEAR_VIBRATIONS_RANGE [AXIS=<x|y|a|b>] [VMIN=<minimal velocity>] [VMAX=<maximal velocity>] [STEP=<steps size of veloctity changes>] `**


## Installation
This extension requires matplotlib. To install it login into your rpi using ssh and excecute the following command:
```~/klippy-env/bin/pip install -v matplotlib```
This will take a couple of minutes to run through and create some load on the rpi. Don't do this while printing. 

Move the file `linear_movement_vibrations.py` (you can find it here in this repo) into the `klippy/extras` directory to install the extension. 

Add the section `[linear_movement_vibrations]` to your `printer.cfg`
Example configuration:
```
[linear_movement_vibrations]
accel_chip: adxl345
x_min: 5
x_max: 175
y_min: 5
y_max: 175
output_directory: /home/pi/klipper_config/linear_vibrations/
```
Make sure the defined output directory is writable. The one in the example configuration shown above will create a folder `linear_vibrations` that can be accessed via the file browser in the web frontend, assuming you are using RatOS. Similar to the `input_shaper folder`, you can find it in the machine tab in the `config` root.  If you are not using RatOS and are unsure which directroy to use `/tmp/` is a save bet. Be aware, that the pngs will not be automatically removed.
`x_min, x_max,y_min, y_max`define a triangle in which the measuremnts will be performend. 
z
If you are using multiple accelerometers, you can also define them as such (untested feature)
```
accel_chip_x: adxl345 rpi
accel_chip_y: adxl345
```
**After these steps, reboot your raspberry, as otherwhise the newly installed pythin libraries are not loaded into the env.**  
