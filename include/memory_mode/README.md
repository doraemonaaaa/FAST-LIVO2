# Native LIVO memory mode

Adapted from `deep_autolabel/nodes/cpp/src/frontend/fast-livo2/memory_mode`.
The native algorithm retains its existing grayscale reference images and
per-frame PCD/RGB export; it does not adopt the migrated photometric pyramid or
its accumulated RGB-map output.

Set the global ROS parameter `memory_mode` to `full` (default), `page`, or
`spatial_sliding`. The last is named `sliding` in the reference implementation;
this native implementation uses the explicitly requested `spatial_sliding` name.
For Scene-Reconstruction, configure it in the runtime YAML:

```yaml
fast_livo2:
  ros_parameters:
    memory_mode: page
```

The implementation lives in `include/memory_mode/{page,spatial_sliding}` and
`src/memory_mode/`; the former standalone `page/` directories were moved, with
no duplicate implementation retained.

| Mode | Matching history | Storage |
|---|---|---|
| `full` | Complete | Resident |
| `page` | Complete | Resident working set plus lossless disk backing |
| `spatial_sliding` | Local spatial window | Resident geometry; lossless image backing; outside history is deleted |

## Spatial sliding

```yaml
fast_livo2:
  ros_parameters:
    memory_mode: spatial_sliding
    sliding_half_extent_m: 50.0
    sliding_trigger_m: 5.0
```

The IMU position centres a world-aligned cube, initially 100 m wide. It
recentres when translation reaches the 5 m trigger. Both distances must be
finite, positive and satisfy trigger < half extent. LIO and every camera share
these world bounds, retaining each voxel that intersects the closed cube,
including boundary cells. Whole retained trees, covariances and observations
are preserved. New voxels wholly outside the current bounds are rejected.

This policy starts after motion/staged initialization and all startup replay
inputs. Pruning runs serially after the complete LIO/VIO/export event: for
image-clock LIVO after the VIO half, and for scan-clock/LiDAR-only input after
each completed event. All cameras' borrowed references are cleared before
deleting any shared visual points or LIO roots. The legacy
`local_map.map_sliding_en` option cannot be combined with this mode.

Removed geometry is deleted rather than paged. Returning to an old area builds
it again. This can change future matches and trajectories; it is not an
equivalent-storage replacement for full/page. Exported scans/images still cover
the whole episode. It is a spatial policy, not a time/observation-count cap:
old observations inside the cube can remain and RSS is not strictly bounded.
Retained grayscale reference images use the same exact-byte backing as `page`:
cold pixels can be reclaimed and faulted back in, while spatial pruning alone
still decides which observations survive. LIO voxel paging stays exclusive to
`page`. This avoids keeping thousands of full-resolution camera frames resident
when several cameras share the spatial window; it does not change matching
geometry, image samples, or window bounds.

The memory CSV appends window activation/centre/half extent, move counts,
deleted LIO roots and visual points, and rejected outside-window insertions.
Initialization logs must have zero window moves/evictions. Large-window runs
with no deletions should retain the full-mode outputs exactly.

`full` retains the original image and voxel residency. `page` changes storage:

- `image_backing.h` writes each processed grayscale image's exact bytes to an
  unlinked file and maps it read-only. Every Feature owns the mapping through a
  shared handle, so a reference remains valid even after the store is destroyed.
  After the frame, current/submap observation images remain working; other
  mappings receive `MADV_DONTNEED`. A later pixel read reloads the same bytes.
  There is no image re-encoding, resizing, or feature deletion.
- `voxel_pager.h` saves complete octree/plane/point/counter state to a temporary
  backing file, replaces cold root pointers with null, and keeps the original
  hash-table keys and iteration order. Every EKF iteration restores the full
  point-query neighborhood before OpenMP. Serial VIO plane lookups also restore
  their roots before dereferencing. Current and preceding LIO-frame accesses
  remain resident. Dirty records are rewritten; stable slots are reused.
- Matching, measurement order, state updates, map point membership,
  sliding-map settings and exported point coordinates follow the same code.
  This is paging, not a smaller matching radius or a new downsampling policy.

Cache files live under the run's `Log/.memory_cache/` (the runner binds Log to
the selected output's `frontend_log`). Image files are unlinked immediately;
mapping ownership controls their lifetime. The voxel backing file is removed
on normal destruction. Records are same-process scratch storage, not portable
maps or restart checkpoints. Disk I/O errors fail the run.

All modes write `memory_working_set.csv` about once per sensor second and once
at exit: process RSS and peak RSS in KiB, total/resident/query voxel counts,
load/eviction/write counters, visual points/observations/reference images, mapped
and working image bytes. The export runner also copies the final pass's CSV to
`raw/memory_working_set.csv`. Initialization subprocesses have separate logs.

Paging does not impose a total RSS cap. Visual feature metadata, hash keys,
sensor queues, current images and the actual query working set remain resident.
The independent batch-bias initializer is unchanged. Publishing the complete
plane map restores all roots for that operation; ordinary pipeline profiles
disable plane publication. Standard pipeline RGB output is already per frame
(`pcd_save.interval: 1`), so it needs no second global RGB spool. Legacy
accumulate-until-exit export remains outside this memory bound.

The offline export runner sets rosbag's publisher queue to 10000 per topic,
with `FLIVO_PLAY_QUEUE_SIZE` as an explicit override. Capacity is not an upfront
allocation; existing trajectory backpressure still controls buffered sensor
time. A small publisher queue can drop IMU messages during paging I/O stalls,
so a valid memory comparison must also check coverage and exported frame counts.

File-backed kernel page cache can remain after process pages are dropped and
can be reclaimed by the OS; it is different from process RSS. Validation must
report measured process peaks and output equality, not only allocated bytes.
