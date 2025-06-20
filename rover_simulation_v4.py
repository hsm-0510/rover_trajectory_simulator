import math
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # Required for 3D plotting

# --- Constants for Earth Model (WGS84) ---
A_SEMI_MAJOR = 6378137.0  # Semi-major axis (equatorial radius) in meters
F_FLATTENING = 1 / 298.257223563  # Flattening
E_SQ = F_FLATTENING * (2 - F_FLATTENING) # Eccentricity squared, derived from flattening

# --- Vector Operations ---
def vec_sub(v1, v2):
    """Subtracts two 3D vectors."""
    return [v1[i] - v2[i] for i in range(3)]

def vec_add(v1, v2):
    """Adds two 3D vectors."""
    return [v1[i] + v2[i] for i in range(3)]

def vec_mul_scalar(v, s):
    """Multiplies a 3D vector by a scalar."""
    return [v[i] * s for i in range(3)]

def vec_magnitude(v):
    """Calculates the Euclidean magnitude of a 3D vector."""
    return math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)

def vec_normalize(v):
    """Normalizes a 3D vector to a unit vector."""
    mag = vec_magnitude(v)
    if mag == 0:
        return [0, 0, 0] # Return zero vector if magnitude is zero
    return [v[i] / mag for i in range(3)]

# --- LLA to ECEF Conversion ---
def lla_to_ecef(lat_deg, lon_deg, alt_m):
    """
    Converts Latitude, Longitude (in degrees), and Altitude (in meters)
    to Earth-Centered, Earth-Fixed (ECEF) Cartesian coordinates (X, Y, Z).
    """
    lat_rad = math.radians(lat_deg)
    lon_rad = math.radians(lon_deg)

    # Radius of curvature in the prime vertical (N)
    N = A_SEMI_MAJOR / math.sqrt(1 - E_SQ * math.sin(lat_rad)**2)

    # Calculate ECEF coordinates
    x = (N + alt_m) * math.cos(lat_rad) * math.cos(lon_rad)
    y = (N + alt_m) * math.cos(lat_rad) * math.sin(lon_rad)
    z = (N * (1 - E_SQ) + alt_m) * math.sin(lat_rad)

    return [x, y, z]

# --- Main Trajectory Generation Function ---
def generate_trajectory(
    start_lat, start_lon, start_alt,
    end_lat, end_lon, end_alt,
    max_speed, max_acceleration, rate_of_acceleration, rate_of_deceleration,
    stationary_duration, sample_time=0.1
):
    """
    Generates a rover's trajectory from a starting LLA point to an ending LLA point
    in ECEF coordinates, considering various motion constraints.
    """
    trajectory_points = []

    # Point of Deceleration
    time_to_decelerate = max_acceleration/rate_of_deceleration
    speed_at_dec = ((max_speed) - (0.5 * (rate_of_deceleration * time_to_decelerate * time_to_decelerate)))
    
    # Convert start and end LLA coordinates to ECEF
    start_ecef = lla_to_ecef(start_lat, start_lon, start_alt)
    end_ecef = lla_to_ecef(end_lat, end_lon, end_alt)

    # Calculate the total linear distance between start and end ECEF points
    displacement_vector_total = vec_sub(end_ecef, start_ecef)
    total_target_distance = vec_magnitude(displacement_vector_total)

    if total_target_distance < 1e-6:
        num_stationary_steps = int(stationary_duration / sample_time)
        for i in range(num_stationary_steps + 1):
            trajectory_points.append([i * sample_time] + start_ecef + [0.0, 0.0])
        return trajectory_points

    direction_unit_vector = vec_normalize(displacement_vector_total)

    current_time = 0.0
    current_position = list(start_ecef)
    current_velocity_magnitude = 0.0
    current_acceleration_magnitude = 0.0
    current_path_distance = 0.0
    current_phase = "stationary"

    max_sim_time = 1000000.0

    while True:
        trajectory_points.append([current_time] + current_position + 
                                 [current_velocity_magnitude, current_acceleration_magnitude])

        if current_path_distance >= total_target_distance - 1e-6:
            if trajectory_points[-1][1:4] != end_ecef:
                if current_velocity_magnitude > 1e-9:
                    time_to_reach_exact_end = (total_target_distance - (current_path_distance - (current_velocity_magnitude * sample_time))) / (current_velocity_magnitude * sample_time) * sample_time
                    current_time = trajectory_points[-1][0] - sample_time + time_to_reach_exact_end
                
                trajectory_points.append([current_time] + end_ecef + [0.0, 0.0])
            break

        dt = sample_time
        next_velocity_magnitude = current_velocity_magnitude
        next_acceleration_magnitude = current_acceleration_magnitude
        
        if current_phase == "stationary":
            if current_time + dt >= stationary_duration:
                current_time = stationary_duration
                current_phase = "jerk_ramp_up"
            else:
                current_time += dt
                continue

        if current_phase == "jerk_ramp_up":
            next_acceleration_magnitude += rate_of_acceleration * dt
            if next_acceleration_magnitude >= max_acceleration:
                next_acceleration_magnitude = max_acceleration
                current_phase = "constant_accel"
            
            avg_accel_for_step = (current_acceleration_magnitude + next_acceleration_magnitude) / 2.0
            next_velocity_magnitude = current_velocity_magnitude + avg_accel_for_step * dt
            
                
            if next_velocity_magnitude >= max_speed:
                next_velocity_magnitude = max_speed
                current_phase = "constant_speed"
                next_acceleration_magnitude = 0.0
            
            if (next_velocity_magnitude >= speed_at_dec) and (next_velocity_magnitude < max_speed):
                next_acceleration_magnitude = current_acceleration_magnitude - (rate_of_deceleration * dt)
                if next_acceleration_magnitude < 0:
                    next_acceleration_magnitude = 0
                next_velocity_magnitude = current_velocity_magnitude + (next_acceleration_magnitude * dt)
                current_phase = "deceleration"
                
        elif current_phase == "constant_accel":
            
            if (next_velocity_magnitude >= speed_at_dec) and (next_velocity_magnitude < max_speed):
                next_acceleration_magnitude = current_acceleration_magnitude - (rate_of_deceleration * dt)
                if next_acceleration_magnitude < 0:
                    next_acceleration_magnitude = 0
                next_velocity_magnitude = current_velocity_magnitude + (next_acceleration_magnitude * dt)
                current_phase = "deceleration"
            
            if next_velocity_magnitude >= max_speed:
                next_velocity_magnitude = max_speed
                current_phase = "constant_speed"
                next_acceleration_magnitude = 0.0
            
            next_acceleration_magnitude = max_acceleration
            next_velocity_magnitude = current_velocity_magnitude + max_acceleration * dt

        elif current_phase == "constant_speed":
            next_velocity_magnitude = max_speed
            next_acceleration_magnitude = 0.0

        elif current_phase == "deceleration":
            next_acceleration_magnitude = current_acceleration_magnitude - (rate_of_deceleration * dt)
            if next_acceleration_magnitude < 0:
                next_acceleration_magnitude = 0
            next_velocity_magnitude = current_velocity_magnitude + (next_acceleration_magnitude * dt)
            
#             avg_accel_for_step = (current_acceleration_magnitude + next_acceleration_magnitude) / 2.0
#             next_velocity_magnitude = current_velocity_magnitude + avg_accel_for_step * dt
            
            if next_velocity_magnitude >= max_speed:
                next_velocity_magnitude = max_speed
                current_phase = "constant_speed"
                next_acceleration_magnitude = 0.0
            
        avg_velocity_for_step = (current_velocity_magnitude + next_velocity_magnitude) / 2.0
        distance_to_move_this_step = avg_velocity_for_step * dt

        projected_current_path_distance = current_path_distance + distance_to_move_this_step

        if projected_current_path_distance >= total_target_distance - 1e-9:
            remaining_dist_for_exact_stop = total_target_distance - current_path_distance
            if avg_velocity_for_step > 1e-9:
                time_fraction_of_step = remaining_dist_for_exact_stop / avg_velocity_for_step
            else:
                time_fraction_of_step = 0.0

            current_time += (time_fraction_of_step * dt)
            current_position = list(end_ecef)
            current_path_distance = total_target_distance
            trajectory_points.append([current_time] + current_position + [0.0, 0.0])
            break
        
        position_change_vector = vec_mul_scalar(direction_unit_vector, distance_to_move_this_step)
        current_position = vec_add(current_position, position_change_vector)
        current_path_distance = projected_current_path_distance
        
        current_time += dt
        current_velocity_magnitude = next_velocity_magnitude
        current_acceleration_magnitude = next_acceleration_magnitude
        
        if current_time > max_sim_time:
            if total_target_distance > 0 and abs(current_path_distance - total_target_distance) > 1e-3:
                trajectory_points.append([current_time] + end_ecef + [0.0, 0.0])
            break

    return trajectory_points

# --- Apply Moving Average Filter ---
def apply_moving_average(data, window_size):
    """Applies a moving average filter to the speed and acceleration data."""
    smoothed_data = []
    for i in range(len(data)):
        if i < window_size - 1:
            # For the initial points, use the average of available points
            start_idx = 0
        else:
            start_idx = i - window_size + 1
        window = data[start_idx:i + 1]
        avg_speed = sum(point[4] for point in window) / len(window)
        avg_accel = sum(point[5] for point in window) / len(window)
        smoothed_data.append([data[i][0]] + data[i][1:4] + [avg_speed, avg_accel])
    return smoothed_data

# --- Write Trajectory Data to CSV File ---
def write_to_csv(filename, data):
    """Writes the trajectory data (time, x, y, z) to a CSV file."""
    try:
        with open(filename, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            for row in data:
                csv_writer.writerow(row[:4]) # Only write time, x, y, z
        print(f"Trajectory successfully saved to '{filename}'.")
    except IOError as e:
        print(f"Error writing to file {filename}: {e}")

# --- Example Usage (Main Execution Block) ---
def main():
    """Defines example parameters, runs the trajectory generation,
    applies smoothing, writes to CSV, and plots the trajectory, speed, and acceleration."""
    print("--- Rover Trajectory Simulation ---")
    
    start_lat = 25.4900
    start_lon = 46.7210
    start_alt = 10.0

    end_lat = 25.4960
    end_lon = 50.721
    end_alt = 18000.0

    max_speed = 1500.0
    max_acceleration = 10
    rate_of_acceleration = 1
    rate_of_deceleration = 0.1
    stationary_duration = 80
    sample_time = 0.1

    output_filename = "simulated_trajectory.csv"

    print(f"\nStarting Point: LLA({start_lat}, {start_lon}, {start_alt})")
    print(f"Ending Point:   LLA({end_lat}, {end_lon}, {end_alt})")
    print(f"Max Speed: {max_speed} m/s, Max Accel: {max_acceleration} m/s^2, Jerk: {rate_of_acceleration} m/s^3")
    print(f"Stationary Duration: {stationary_duration} s, Sample Time: {sample_time} s")

    print("\nGenerating trajectory...")
    trajectory_data = generate_trajectory(
        start_lat, start_lon, start_alt,
        end_lat, end_lon, end_alt,
        max_speed, max_acceleration, rate_of_acceleration, rate_of_deceleration,
        stationary_duration, sample_time
    )

    if trajectory_data:
        # Apply moving average filter with a window size of 5
        window_size = 5
        smoothed_trajectory = apply_moving_average(trajectory_data, window_size)

        times = [point[0] for point in smoothed_trajectory[:-1]]
        xs = [point[1] for point in smoothed_trajectory[:-1]]
        ys = [point[2] for point in smoothed_trajectory[:-1]]
        zs = [point[3] for point in smoothed_trajectory[:-1]]
        speeds = [point[4] for point in smoothed_trajectory[:-1]]
        accelerations = [point[5] for point in smoothed_trajectory[:-1]]

        write_to_csv(output_filename, smoothed_trajectory[:-1])
        print(f"\nTotal trajectory points generated: {len(smoothed_trajectory)}")
        print(f"First trajectory point: {smoothed_trajectory[0]}")
        print(f"Last trajectory point:  {smoothed_trajectory[-1]}")
        print(f"Simulation completed in approximately {smoothed_trajectory[-1][0]:.2f} seconds.")

        plt.style.use('seaborn-v0_8-darkgrid')

        fig_3d = plt.figure(figsize=(10, 8))
        ax_3d = fig_3d.add_subplot(111, projection='3d')
        ax_3d.plot(xs, ys, zs, label='Rover Trajectory', color='blue', linewidth=2)
        ax_3d.scatter(xs[0], ys[0], zs[0], color='green', s=100, label='Start Point', marker='o')
        ax_3d.scatter(xs[-1], ys[-1], zs[-1], color='red', s=100, label='End Point', marker='X')
        ax_3d.set_xlabel('ECEF X (m)')
        ax_3d.set_ylabel('ECEF Y (m)')
        ax_3d.set_zlabel('ECEF Z (m)')
        ax_3d.set_title('Rover Trajectory in ECEF Coordinates')
        ax_3d.legend()
        ax_3d.grid(True)

        plt.figure(figsize=(10, 6))
        plt.plot(times, speeds, color='purple', linewidth=2)
        plt.xlabel('Time (s)')
        plt.ylabel('Net Speed (m/s)')
        plt.title('Rover Net Speed Over Time')
        plt.grid(True)
        plt.axhline(y=max_speed, color='red', linestyle='--', label='Max Speed Limit')
        plt.legend()

        plt.figure(figsize=(10, 6))
        plt.plot(times, accelerations, color='orange', linewidth=2)
        plt.xlabel('Time (s)')
        plt.ylabel('Net Acceleration (m/s²)')
        plt.title('Rover Net Acceleration Over Time')
        plt.grid(True)
        plt.axhline(y=max_acceleration, color='red', linestyle='--', label='Max Acceleration Limit')
        plt.legend()

        plt.tight_layout()
        plt.show()

    else:
        print("No trajectory points were generated. Check input parameters.")

if __name__ == "__main__":
    main()
