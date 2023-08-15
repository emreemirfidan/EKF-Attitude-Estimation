# EKF Quaternion Based Attitude Estimation
This is a quaternion-based Extended Kalman filter for real-time estimation of the orientation of a UAV.

I tested with the STEVAL-STLCS02V1 (SensorTile) evolution board. I used the lsm303agr accelerometer-magnetometer and lsm6dsm gyroscope sensors on the board. I tested it with the STM32L476JG microcontroller that this board is equipped with.

## How To Use
```c
// In The Main Function:
ekf_t ekf;
float euler_ekf[3];
EKF_init(&ekf, mag.x, mag.y, mag.z, 0.1, 1, 100);
```

```c
// In The While Loop
// Run this function at every delta_time interval (10ms)
EKF_update(&ekf, euler_ekf, accel.x, accel.y, accel.z, gyro.x, gyro.y, gyro.z, mag.x, mag.y, mag.z, delta_time_second);

// Convert radians to degrees
roll_deg = euler_ekf[0] * 180.0f / 3.14159265f;
pitch_deg = euler_ekf[1] * 180.0f / 3.14159265f;
yaw_deg = euler_ekf[2] * 180.0f / 3.14159265f;
```
