# Purpose

The purpose of this code is to simulate a rover trajectory from a starting geographic location to an ending location, incorporating a stationary phase, increasing acceleration phase, constant acceleration phase, decreasing acceleration phase, and constant speed phase. The simulation generates a 3D trajectory, calculates and plots net velocity and net acceleration, and saves the trajectory data to a CSV file through a moving average filter to smooth the transitions between the phases in order to avoid sudden unwanted jerks. This trajectory generator is directly applicable for the gps-sdr-sim GPS Simulator by osqzss, and can be performed through these steps:

# Generate The Simulation

1. In the code find the def main() and adjust the start_lat, start_lon, start_alt, end_lat, end_lon, end_alt, max_speed, max_acceleration, rate_of_acceleration, rate_of_deceleration, stationary_duration, and sample_time accoring to your trajectory requirements.
2. Run the python code rover_simulation_v4.py and look for the generated csv file 'simulated_trajectory.csv'
3. Copy the csv file at the installation directory of gps-sdr-sim and run the following command:

```
./gps-sdr-sim -e brdc3610.24n -u simulated_trajectory.csv -b 8
```

Note: the -e brdc3160.24 can change according to the ephemeris file being used, the ephemeris files can be downloaded through link provided at gps-sdr-sim repository.

4. The gpssim.bin file would be created and ready for trajectory simulation.

See following link for gps-sdr-sim installation guide: https://github.com/osqzss/gps-sdr-sim

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

4. Generating position data by integrating a smoothed speed profile and applying a linear height offset.

5. Saving the resulting x, y, z coordinates along with time to a CSV file.

6. Plotting the 3D trajectory and analyzing net velocity and acceleration, with both raw and smoothed velocity profiles displayed for comparison. The code is designed to be modular, allowing adjustments to parameters like window_size and ramp_duration to fine-tune the smoothing and transition effects.


# Example Parameters
```
    start_lat = 25.4900
    start_lon = 46.7210
    start_alt = 10.0

    end_lat = 25.4960
    end_lon = 50.721
    end_alt = 18000.0

    max_speed = 1500.0
    max_acceleration = 10
    rate_of_acceleration = 1
    rate_of_deceleration = 1
    stationary_duration = 80
    sample_time = 0.1

    window_size = 350
```
# 3D Plot

![image](https://github.com/user-attachments/assets/e72df3d5-43a4-4149-959f-26ef1b6757d1)


# Net Velocity Plot

![image](https://github.com/user-attachments/assets/45c9ecfd-e045-4b11-9db9-9f7cd242d5af)


# Net Acceleration Plot

![image](https://github.com/user-attachments/assets/43f3d782-cb11-46b1-b31a-8380034da694)

---

# Future Works

1. Problem: Sudden Jerks appear when selected rate of acceleration or deceleration do not match with final velocity, which leads to premature halt at max speed and suddent jerk.
![image](https://github.com/user-attachments/assets/39efbea4-3779-4b7f-9147-0822b9e3ebd1)
2. Solution: Compute the rate of deceleration that achieves longest deceleration time without jerk. If the computed rate of deceleration is more than user's, then implement the computed rate of deceleration. However, if the computed rate of deceleration is less than user's, then implement the user's rate of deceleration instead.
3. Problem: Increase number of waypoints for trajectory simulation with non-jerk directional changes
4. Solution: Use of polynomial equation-based path estimation for non-jerk movement through the user-defined waypoints
