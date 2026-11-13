# Initialization modules

Headers live in [`include/initialization/`](../../include/initialization/).
Source files live in this directory. Online IMU propagation and ordinary LIO/VIO
updates remain in their existing modules.

The default uses independent batch registration again. The experimental native
LIO reuse has been reverted; the historical temporary LIVO clock window remains.

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
LiDAR-only, frozen once and excludes query scans. Clock scoring reports activity
through the existing replay audit so the offline runner does not abort it as
an idle EOF window.

The external runner `our_method/run_fast_livo2_export.bash` launches fresh
processes. `our_method/pipeline/staged_initialization.py` transfers means and
clock decisions; `initialization_pass.py` binds input/config/runtime identities.
Configuration and limitations:
[`our_method/README.md`](../../../../our_method/README.md#moving-start-deeprobotics-bags).

`imu_support.h` is shared by initialization and normal processing. Batch fitting
uses fully covered, continuous IMU/LiDAR segments; unsupported scans are excluded
before independent registration. Motion endpoints are bracketed by real IMU,
and rewind never extends boundary measurements. Normal scheduling rejects entire
unsupported scans before any visual slicing; both replay passes select identical
IMU brackets from sensor event times rather than filter initialization state.

## Parallel initialization tasks

`imu/initialization_threads` controls initialization tasks (default 4 in the
DeepRobotics method config; 1 executes the same schedule serially). Normal
per-frame LIO/VIO threading remains independent of this setting. OpenMP runtime
limits can cap the requested worker count.

- Stage one prepares each target KD-tree/normal geometry independently. Pairs
  with a valid previous-pass guess register concurrently. The first pass and
  pairs without that guess retain the current-pass preceding-pair fallback,
  executed in input order after independent pairs finish.
- A pass freezes its registration results before independent per-scan LiDAR
  deskew. The task barrier completes all clouds before the next pass.
- Global poses/longest continuous segments accumulate in input order. Ceres
  and Eigen remain single-threaded; bias means and prior covariance semantics
  are unchanged.
- Stage-two early motion fitting and the three-second temporary LIVO sequence
  remain chronological. Offset reference registration uses the same pair
  scheduler, followed by independent scan deskew, final pair refits and frozen
  query-map construction.
- Offset candidates remain sequential. Within a candidate, query IMU deskew,
  registration and scoring run independently. Score reductions and CSV output
  preserve input order. Each query owns its frozen map; workers do not query a
  shared mutable KD-tree.

`parallel.hpp` keeps worker writes in distinct preallocated output slots and
rethrows worker exceptions in input order on the caller. Progress callbacks
also execute only on the caller thread, including the offline replay heartbeat.

`batch_initialization_timing.json` separates input preparation, independent
LiDAR trajectory construction and the single-thread batch solve.
`imu_time_offset_timing.json` reports reference construction, frozen query maps
and candidate scoring. These sidecars do not change algorithm diagnostics or
selection policy. The original sequential fallback schedule and one/four-worker
execution were validated with bitwise candidate-score/CSV comparisons before
removing the regression sources and targets. Results are recorded in the
[initialization parallelization experiment](../../../../output/experiments/initialization_parallel_20261113_01/README.md).
