//! Updraft helicity diagnostic.

use rayon::prelude::*;

use crate::compute::ComputeOpts;
use crate::error::WrfResult;
use crate::file::WrfFile;

/// Linearly interpolate value at height `z` between two levels.
/// `z0`/`z1` are heights at the two levels, `v0`/`v1` are the field values.
#[inline]
fn lerp_at(z: f64, z0: f64, z1: f64, v0: f64, v1: f64) -> f64 {
    if (z1 - z0).abs() < 1e-6 {
        0.5 * (v0 + v1)
    } else {
        let frac = (z - z0) / (z1 - z0);
        v0 + frac * (v1 - v0)
    }
}

/// Updraft helicity (m^2/s^2). `[ny, nx]`
///
/// UH = integral from z_bot to z_top of (w * zeta_z) dz
/// Default layer: 2-5 km AGL.
///
/// The integration uses trapezoidal quadrature.  When a model layer only
/// partially overlaps the integration bounds, w and vorticity are linearly
/// interpolated to the exact boundary height so that the partial-layer
/// contribution is accurate.  Only positive vertical velocity is counted
/// (updraft helicity), matching WRF's internal calculation.
pub fn compute_uhel(f: &WrfFile, t: usize, opts: &ComputeOpts) -> WrfResult<Vec<f64>> {
    let w = f.w_destag(t)?;
    let u = f.u_destag(t)?;
    let v = f.v_destag(t)?;
    let h_agl = f.height_agl(t)?;

    let nx = f.nx;
    let ny = f.ny;
    let nz = f.nz;
    let nxy = nx * ny;
    let dx = f.dx;
    let dy = f.dy;

    // Integration layer bounds (default 2-5 km AGL), configurable
    let z_bot = opts.bottom_m.unwrap_or(2000.0);
    let z_top = opts.top_m.unwrap_or(5000.0);

    // Compute relative vorticity at each level.
    // vorticity(u, v, ...) computes dv/dx - du/dy.
    let mut vort_3d = vec![0.0f64; nz * nxy];
    vort_3d.par_chunks_mut(nxy).enumerate().for_each(|(k, plane)| {
        let u_plane = &u[k * nxy..(k + 1) * nxy];
        let v_plane = &v[k * nxy..(k + 1) * nxy];
        let vort = crate::met::dynamics::vorticity(u_plane, v_plane, nx, ny, dx, dy);
        plane.copy_from_slice(&vort);
    });

    // Integrate w * vort over the layer for each column using trapezoidal rule.
    // At partial-layer boundaries, interpolate w and vort to the exact bound.
    let mut uhel = vec![0.0f64; nxy];
    uhel.par_iter_mut().enumerate().for_each(|(ij, uh_val)| {
        let mut integral = 0.0f64;
        for k in 0..nz - 1 {
            let idx0 = k * nxy + ij;
            let idx1 = (k + 1) * nxy + ij;

            let h0 = h_agl[idx0];
            let h1 = h_agl[idx1];

            // Skip layers entirely outside the integration bounds
            if h1 <= z_bot || h0 >= z_top {
                continue;
            }

            // Clamp to integration bounds
            let z_lo = h0.max(z_bot);
            let z_hi = h1.min(z_top);
            let dz = z_hi - z_lo;
            if dz <= 0.0 {
                continue;
            }

            // Interpolate w and vort to the clamped endpoints
            let w0_raw = w[idx0];
            let w1_raw = w[idx1];
            let v0_raw = vort_3d[idx0];
            let v1_raw = vort_3d[idx1];

            let w_lo = if z_lo > h0 {
                lerp_at(z_lo, h0, h1, w0_raw, w1_raw)
            } else {
                w0_raw
            };
            let w_hi = if z_hi < h1 {
                lerp_at(z_hi, h0, h1, w0_raw, w1_raw)
            } else {
                w1_raw
            };
            let vort_lo = if z_lo > h0 {
                lerp_at(z_lo, h0, h1, v0_raw, v1_raw)
            } else {
                v0_raw
            };
            let vort_hi = if z_hi < h1 {
                lerp_at(z_hi, h0, h1, v0_raw, v1_raw)
            } else {
                v1_raw
            };

            // Only count where w > 0 (updraft)
            let w_avg = 0.5 * (w_lo + w_hi);
            if w_avg > 0.0 {
                let vort_avg = 0.5 * (vort_lo + vort_hi);
                integral += w_avg * vort_avg * dz;
            }
        }
        *uh_val = integral;
    });

    Ok(uhel)
}
