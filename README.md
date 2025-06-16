# Purpose

The purpose of this code is to simulate a rover trajectory from a starting geographic location to an ending location, incorporating a stationary phase, acceleration, constant speed, and deceleration phases, while adding a parabolic height profile. The simulation generates a 3D trajectory, calculates and plots net velocity and net acceleration, and saves the trajectory data to a CSV file for further analysis or visualization.

# Design Method

The code employs a physics-based approach to model the rover's motion:

1. Coordinate Transformation: Utilizes the pyproj.Transformer to convert geodetic coordinates (latitude, longitude, altitude) to Earth-Centered Earth-Fixed (ECEF) coordinates for 3D spatial representation.

2. Motion Profile: Defines a motion profile with a stationary duration, followed by an acceleration phase, a constant speed phase, and a deceleration phase. The distances and times are scaled to match the total distance between start and end points.

3. Smoothing: Applies a moving average filter to smooth the velocity profile, ensuring a more realistic transition between motion phases, with a customizable window size and ramp duration for gradual speed changes.

4. Integration: Uses numerical integration (trapezoidal rule) of the smoothed speed profile to compute distances, which are then used to interpolate positions between start and end points.

5. Visualization: Employs matplotlib to create 3D trajectory plots and 2D plots of net velocity and acceleration over time.

6. Data Storage: Saves the simulated trajectory (time, x, y, z coordinates) to a CSV file without headers for compatibility with external tools.

# Description

This Python script simulates a rover's trajectory by:

1. Defining input parameters such as start/end coordinates, maximum height, speed, acceleration, and time durations.

2. Converting geographic coordinates to ECEF using pyproj.

3. Calculating a scaled motion profile to fit the total distance, including a stationary period and smoothed acceleration/deceleration phases.

4. Generating position data by integrating a smoothed speed profile and applying a parabolic height offset.

5. Saving the resulting x, y, z coordinates along with time to a CSV file.

6. Plotting the 3D trajectory and analyzing net velocity and acceleration, with both raw and smoothed velocity profiles displayed for comparison. The code is designed to be modular, allowing adjustments to parameters like window_size and ramp_duration to fine-tune the smoothing and transition effects.
