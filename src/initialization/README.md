# Initialization modules

Headers live in [`include/initialization/`](../../include/initialization/).
Source files live in this directory. Online IMU propagation and ordinary LIO/VIO
updates remain in their existing modules.

| Responsibility | Files |
|---|---|
| Independent LiDAR trajectory, three registration/deskew passes | `batch_lidar_trajectory.cpp` |
| Joint preintegrated IMU/pose/velocity/gravity/bias batch optimizer | `batch_imu_initializer.hpp`, `batch_imu_initializer.cpp` |
| ROS1 bag reader using the native preprocessor; bias-only handoff | `batch_main.cpp`, `input.h` |
| Motion velocity/gravity fitting with batch bias seeds | `motion_initializer.h` |
| Temporary post-motion state window and clock decision handoff | `staged_initialization.h`, `staged_initialization.cpp` |
| Frozen LiDAR reference, candidate deskew, balanced score and holdout policy | `imu_time_offset/` |
| Soft gravity-norm model, correlated state/Joseph covariance update | `gravity_norm_constraint.h`, `gravity_norm_constraint.cpp` |
| Offline checkpoint lifecycle, full-state serialization and reread verification | `startup_replay.*`, `initialization_checkpoint.h` |
| Backward integration and IMU propagation context | `imu_rewind.cpp` |
| Reference diagnostics (no acceptance gates) | `initialization_diagnostics.h` |

```text
independent LiDAR/IMU batch optimization
    -> bias_g and bias_a MEANS ONLY; original ESIKF covariance
    -> reread and fit early velocity/gravity
    -> temporary LIVO for 3 seconds, biases continue updating
    -> frozen clock scoring (optional)
    -> if clock changed: independent batch bias calibration on corrected times
    -> refit early motion using batch means
    -> rewind/checkpoint, fresh maps and formal replay
    -> normal LIVO, with optional soft gravity-norm model
```

The standalone stage-one executable creates no `LIVMapper` or ESIKF. Its batch
velocity/gravity/covariance are not transferred. Stage two uses the ordinary
shared LIVO state. `LIVMapper` calls `StartupReplay`, which owns
`StagedInitialization`; complete raw scans are retained only for the short clock
window, irrespective of the visual event scheduler. Clock reference geometry is
LiDAR-only, frozen once and excludes query scans.

The external runner `our_method/run_fast_livo2_export.bash` launches fresh
processes. `our_method/pipeline/staged_initialization.py` transfers means and
clock decisions; `initialization_pass.py` binds input/config/runtime identities.
Configuration and limitations:
[`our_method/README.md`](../../../../our_method/README.md#moving-start-deeprobotics-bags).
