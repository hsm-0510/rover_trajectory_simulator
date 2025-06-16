import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyproj import Transformer

# =======================
# INPUT PARAMETERS
# =======================
start_lat = 25.4900              # degrees
start_lon = 46.7210              # degrees
start_alt = 10                    # meters

end_lat = 25.5230                # degrees
end_lon = 50.4610                # degrees
end_alt = 10                      # meters

max_height = 18000              # meters above the highest of start/end alt
max_speed = 1500                 # meters per second
max_acceleration = 40            # meters per second squared
duration_at_max_speed = 60      # seconds
stationary_duration = 80         # seconds
time_step = 0.1                  # seconds
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
    max_height, max_speed, max_accel, duration_at_max_speed,
    stationary_duration, time_step, output_csv
):
    # Convert start and end positions to ECEF
    start_ecef = geodetic_to_ecef(start_lat, start_lon, start_alt)
    end_ecef = geodetic_to_ecef(end_lat, end_lon, end_alt)
    total_distance = np.linalg.norm(end_ecef - start_ecef)

    # Motion profile setup
    t_accel = max_speed / max_accel
    d_accel = 0.5 * max_accel * t_accel**2
    d_const = max_speed * duration_at_max_speed
    d_total = 2 * d_accel + d_const

    # Scale factor to adjust actual distance to match required
    scale = total_distance / d_total

    # Adjust times and distances by scale
    t_accel *= np.sqrt(scale)
    max_speed *= np.sqrt(scale)
    t_const = duration_at_max_speed * scale
    d_accel = 0.5 * max_accel * t_accel**2
    d_const = max_speed * t_const

    # Total simulation time
    t_total = stationary_duration + 2 * t_accel + t_const
    times = np.arange(0, t_total + time_step, time_step)

    positions = []

    for t in times:
        if t < stationary_duration:
            pos = start_ecef
        else:
            t_rel = t - stationary_duration

            if t_rel < t_accel:
                d = 0.5 * max_accel * t_rel**2
            elif t_rel < t_accel + t_const:
                d = d_accel + max_speed * (t_rel - t_accel)
            else:
                t_dec = t_rel - (t_accel + t_const)
                d = d_accel + d_const + max_speed * t_dec - 0.5 * max_accel * t_dec**2

            # Cap distance at total
            fraction = min(d / (d_accel + d_const + d_accel), 1.0)
            pos = (1 - fraction) * start_ecef + fraction * end_ecef

            # Add height offset (parabolic bump)
            h_frac = 4 * fraction * (1 - fraction)
            height_offset = (max_height - max(start_alt, end_alt)) * h_frac
            up_dir = pos / np.linalg.norm(pos)
            pos += up_dir * height_offset

        positions.append(pos)

    positions = np.array(positions)

    # Save to CSV without header
    df = pd.DataFrame({
        'time': times,
        'x': positions[:, 0],
        'y': positions[:, 1],
        'z': positions[:, 2]
    })
    df.to_csv(output_csv, index=False, header=False)

    # Plot trajectory
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], label='Simulated Trajectory')
    ax.set_title("3D Rocket Trajectory")
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.set_zlabel("Z (m)")
    ax.legend()
    plt.tight_layout()
    plt.show()
    
        # Calculate Net Velocity (magnitude of first derivative)
    velocities = np.gradient(positions, time_step, axis=0)
    speed = np.linalg.norm(velocities, axis=1)

    # Calculate Net Acceleration (magnitude of second derivative)
    accelerations = np.gradient(velocities, time_step, axis=0)
    accel_mag = np.linalg.norm(accelerations, axis=1)

    # Plot Net Velocity
    plt.figure(figsize=(10, 4))
    plt.plot(times, speed, label="Net Velocity", color='blue')
    plt.title("Net Velocity vs Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Speed (m/s)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Plot Net Acceleration
    plt.figure(figsize=(10, 4))
    plt.plot(times, accel_mag, label="Net Acceleration", color='red')
    plt.title("Net Acceleration vs Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Acceleration (m/s²)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    return df


# =======================
# RUN SIMULATION
# =======================
simulate_trajectory(
    start_lat, start_lon, start_alt,
    end_lat, end_lon, end_alt,
    max_height, max_speed, max_acceleration,
    duration_at_max_speed, stationary_duration,
    time_step, output_csv
)
