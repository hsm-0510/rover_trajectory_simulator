import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyproj import Transformer

# =======================
# INPUT PARAMETERS
# =======================
start_lat = 25.4900       # degrees
start_lon = 46.7210       # degrees
start_alt = 10            # meters

end_lat = 25.4960         # degrees
end_lon = 50.0000         # degrees
end_alt = 18000           # meters

max_height = end_alt      # meters above the highest of start/end alt
max_speed = 1500        # meters per second  <-- This is the ceiling
max_acceleration = 40   # meters per second squared
stationary_duration = 80   # seconds
time_step = 0.1         # seconds
output_csv = "simulated_trajectory.csv"  # output filename

# =======================
# FUNCTION DEFINITIONS
# =======================


def geodetic_to_ecef(lat, lon, alt):
    transformer = Transformer.from_crs("epsg:4979", "epsg:4978", always_xy=True)
    x, y, z = transformer.transform(lon, lat, alt)
    return np.array([x, y, z])

def simulate_trajectory(
    start_lat, start_lon, start_alt,
    end_lat, end_lon, end_alt,
    max_height, max_speed, max_accel,
    stationary_duration, time_step, output_csv
):
    # Convert start and end positions to ECEF
    start_ecef = geodetic_to_ecef(start_lat, start_lon, start_alt)
    end_ecef = geodetic_to_ecef(end_lat, end_lon, end_alt)
    total_distance_ecef = np.linalg.norm(end_ecef - start_ecef)

    # --- Calculate motion profile durations and distances based on actual max_speed ---

    # Time to reach max_speed with max_accel
    t_accel_to_max_speed = max_speed / max_accel
    # Distance covered during acceleration to max_speed
    d_accel_to_max_speed = 0.5 * max_accel * t_accel_to_max_speed**2

    # Check if the total_distance_ecef is even enough to reach max speed
    if total_distance_ecef <= d_accel_to_max_speed: # If total_distance is less than accel distance to max speed
        # This implies we might not even reach max_speed
        t_accel_actual = np.sqrt(2 * total_distance_ecef / max_accel)
        d_accel_actual = total_distance_ecef
        t_const = 0 # No constant speed phase
        max_speed_actual = max_accel * t_accel_actual
    else:
        # We accelerate to max_speed and then hold it for the remaining distance
        t_accel_actual = t_accel_to_max_speed
        d_accel_actual = d_accel_to_max_speed
        d_remaining = total_distance_ecef - d_accel_actual
        t_const = d_remaining / max_speed if max_speed > 0 else 0
        max_speed_actual = max_speed # We reached the intended max_speed


    # Total simulation time (stationary + acceleration + constant speed)
    t_total = stationary_duration + t_accel_actual + t_const
    times = np.arange(0, t_total + time_step, time_step)
    
    # Moving Average Window-Size
    window_size = 250
    ramp_duration = 80 # Duration for initial gentle ramp-up

    positions = []
    # Pre-compute a smoothed speed profile with gradual ramp-up
    speed_profile = np.zeros_like(times)
    for i, t in enumerate(times):
        t_rel = t - stationary_duration # Time relative to the end of stationary phase

        if t_rel < 0:
            speed_profile[i] = 0 # Still in stationary phase
        elif t_rel < ramp_duration:
            # Smooth ramp-up using a quadratic function for a gentler start
            current_speed = max_accel * t_rel * (t_rel / ramp_duration)**2
            speed_profile[i] = min(current_speed, max_speed_actual)
        elif t_rel < t_accel_actual: # Changed from t_accel_actual + ramp_duration
            # Continue acceleration phase (linear after ramp-up)
            speed_at_ramp_end = max_accel * ramp_duration * (ramp_duration / ramp_duration)**2 # This is max_accel * ramp_duration
            current_speed = speed_at_ramp_end + max_accel * (t_rel - ramp_duration)
            speed_profile[i] = min(current_speed, max_speed_actual)
        else:
            # Maintain max_speed_actual after the acceleration phase
            speed_profile[i] = max_speed_actual

    smoothed_speed_for_pos_calc = np.convolve(speed_profile, np.ones(window_size)/window_size, mode='same')

    # Ensure the smoothed speed doesn't exceed the intended max_speed
    smoothed_speed_for_pos_calc = np.minimum(smoothed_speed_for_pos_calc, max_speed_actual)

    for i, t in enumerate(times):
        if t < stationary_duration:
            pos = start_ecef
        else:
            # Use cumulative integral of smoothed speed for distance
            idx_motion_start = np.where(times >= stationary_duration)[0][0]
            # Make sure we don't try to integrate an empty slice
            if i + 1 <= idx_motion_start:
                 d = 0.0
            else:
                 d = np.trapezoid(smoothed_speed_for_pos_calc[idx_motion_start : i + 1], times[idx_motion_start : i + 1])


            # Cap distance at total_distance_ecef
            fraction = min(d / total_distance_ecef, 1.0)
            pos = (1 - fraction) * start_ecef + fraction * end_ecef

            # Add height offset (parabolic bump)
            if fraction < 0.5:
                h_frac = 2 * fraction  # Linear increase to max height at midpoint
            else:
                h_frac = 2 * (1 - fraction)  # Linear decrease from max height after midpoint
            
            height_offset = (max_height - max(start_alt, end_alt)) * h_frac
            up_dir = pos / np.linalg.norm(pos)
            pos += up_dir * height_offset

        positions.append(pos)

    positions = np.array(positions)

    # --- Calculate Net Velocity and Acceleration BEFORE Clipping ---
    velocities = np.gradient(positions, time_step, axis=0)
    speed = np.linalg.norm(velocities, axis=1)

    # Apply moving average filter for smoothing for the *output plot*
    smoothed_speed_output = np.convolve(speed, np.ones(window_size)/window_size, mode='same')

    # Calculate Net Acceleration (magnitude of second derivative)
    accelerations = np.gradient(velocities, time_step, axis=0)
    accel_mag = np.linalg.norm(accelerations, axis=1)

    # --- Clipping Logic to remove deceleration tail ---
    
    clip_index = len(times) - 1 # Default to no clipping

    # Find the first index from the end where smoothed speed is below max_speed_actual * 0.99
    # (i.e., where it starts to drop from the plateau)
    decel_start_idx = len(times) - 1 # Initialize to the very end
    
    # Iterate backwards from the end of the steady state (or end of array)
    # to find where the speed *last* consistently maintained the max speed.
    # Start checking from after the acceleration phase, if applicable.
    
    # Find the index after which we expect either constant speed or the end of trajectory
    steady_state_check_start_idx = np.where(times >= stationary_duration + t_accel_actual)[0]
    if len(steady_state_check_start_idx) > 0:
        steady_state_check_start_idx = steady_state_check_start_idx[0]
    else:
        steady_state_check_start_idx = 0 # If no acceleration phase, start from beginning

    # Find the last point that is still at max speed (within a tolerance)
    last_max_speed_idx = -1
    for j in range(len(smoothed_speed_output) - 1, steady_state_check_start_idx, -1):
        if smoothed_speed_output[j] >= max_speed_actual * 0.99: # 99% of max speed
            last_max_speed_idx = j
            break
    
    if last_max_speed_idx != -1:
        # We found a point where it's still at max speed.
        # Clip at this point or a few steps after it if the drop is very sudden.
        # Add a small buffer after the last max speed point, but not too much.
        clip_index = min(last_max_speed_idx + int(5/time_step), len(times) - 1) # Add up to 5 seconds buffer

    # General fallback: if no plateau, or if the clipping is too aggressive
    # ensure we don't clip the entire motion.
    if clip_index < stationary_duration + int(ramp_duration / time_step) + 2:
        # If the automatic clipping is too aggressive, just trim the last 10 seconds of data,
        # or less if the simulation is shorter. This is a failsafe.
        clip_index = max(0, len(times) - int(10/time_step) -1)
        # Ensure it doesn't go below the stationary phase
        clip_index = max(clip_index, np.where(times >= stationary_duration)[0][0] + 5) # Ensure at least 5s of motion

    # Ensure clip_index is valid and doesn't remove the entire simulation
    if clip_index >= len(times) -1: # If no significant drop, use the full array
        clip_index = len(times) -1


    times_clipped = times[:clip_index + 1]
    positions_clipped = positions[:clip_index + 1]
    speed_clipped = speed[:clip_index + 1]
    smoothed_speed_output_clipped = smoothed_speed_output[:clip_index + 1]
    accel_mag_clipped = accel_mag[:clip_index + 1]


    # Save to CSV (without header) - use clipped data
    df_clipped = pd.DataFrame({
        'time': times_clipped,
        'x': positions_clipped[:, 0],
        'y': positions_clipped[:, 1],
        'z': positions_clipped[:, 2]
    })
    df_clipped.to_csv(output_csv, index=False, header=False)

    # Plot trajectory - use clipped data
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(positions_clipped[:, 0], positions_clipped[:, 1], positions_clipped[:, 2], label='Simulated Trajectory')
    ax.set_title("3D Rocket Trajectory")
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.set_zlabel("Z (m)")
    ax.legend()
    plt.tight_layout()
    plt.show()
    
    # Plot Net Velocity - use clipped data
    plt.figure(figsize=(10, 4))
    plt.plot(times_clipped, speed_clipped, label="Net Velocity (Raw)", color='blue')
    #plt.plot(times_clipped, smoothed_speed_output_clipped, label="Net Velocity (Smoothed)", color='red')
    plt.axhline(y=max_speed, color='green', linestyle='--', label=f'Desired Max Speed ({max_speed} m/s)')
    plt.title("Net Velocity vs Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Speed (m/s)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Plot Net Acceleration - use clipped data
    plt.figure(figsize=(10, 4))
    plt.plot(times_clipped, accel_mag_clipped, label="Net Acceleration", color='red')
    plt.title("Net Acceleration vs Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Acceleration (m/s²)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    return df_clipped

# =======================
# RUN SIMULATION
# =======================
simulate_trajectory(
    start_lat, start_lon, start_alt,
    end_lat, end_lon, end_alt*1.273757899,
    max_height*1.273757899, max_speed, max_acceleration/2.16, stationary_duration,
    time_step, output_csv
)
