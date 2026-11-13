# Photometric patch model

`VisualCamera` uses the same sampling and residual construction for every configured camera. The shared ESIKF and per-camera exposure states remain the update mechanism.

## Additive brightness changes

For each patch of N pixels, eliminate an unconstrained additive brightness offset analytically:

```
r_i = w [ a_cur (I_cur,i - mean(I_cur)) - a_ref (I_ref,i - mean(I_ref)) ]
H_i = H_raw,i - mean(H_raw)
```

`a_cur` and `a_ref` are the existing inverse exposure gains. The noise weight `w` is frozen once per frame, as before. Removing the patch mean retains texture contrast and the exposure-gain derivative; it does not normalize away contrast. No extra filter state or persistent patch offset is introduced.

Both residuals and Jacobians must be centered. This is the projection obtained by eliminating a nuisance additive offset from least squares. The constant brightness mode carries no pose information, and a constant patch consequently contributes zero information. Patch admission and initial reference selection use the same mean-centered squared difference, so a uniform brightness shift is not rejected before the update. NCC reference lifecycle selection remains mean invariant.

The optional inverse-compositional update also centers its residual and Jacobian. It retains its existing six-column pose model without estimating exposure gains.

## Coarse sampling

At level `l` and search level `s`, sample locations are:

```
q(x,y) = projected_center + (x - half_patch, y - half_patch) * 2^(l+s)
```

Each location uses the **adjacent four original-image pixels** for bilinear interpolation, including reference affine warping. The coarse stride separates sample locations; it does not set the distance between the four interpolation neighbours. Admission at level zero also samples the current patch at its selected search stride, matching the warped reference patch. Sampling respects the image row stride, including non-contiguous OpenCV ROIs.

The forward Jacobian uses the analytic derivatives of that same bilinear interpolant, in original-image pixel units. It does not divide them again by the coarse stride. The IC reference cache uses the same interpolant and level/search spacing. The conservative border guard is retained.

This remains sparse coarse sampling of the original grayscale image. It does not introduce a Gaussian image pyramid or an anti-aliasing filter.

## Validation history

Validation covered the actual forward and optional IC builders, coarse/search levels, ROI stride, unclipped +/-30 brightness shifts, pose/exposure finite differences and flat patches, plus shared-state/camera-specific Jacobians and dispatch. The regression sources and build targets have since been removed; the [initialization parallelization experiment](../../output/experiments/initialization_parallel_20261113_01/README.md) retains the validation records.

Uniform additive shifts without clipping are invariant under this model. Saturation, moving shadows within a patch, occlusion and incorrect geometry are not removed by subtracting a single mean.
