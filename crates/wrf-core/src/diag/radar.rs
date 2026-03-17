//! Radar diagnostic variables: dbz, maxdbz
//!
//! Simulated reflectivity from hydrometeor mixing ratios following
//! WRF's `module_diag_functions.F` (`REFL_10CM`).
//!
//! The equivalent reflectivity factor (mm^6 m^-3) for each species is:
//!
//!   Z_e = C * (rho_air * q)^1.75
//!
//! where C = factor / (pi^1.75 * rho_x^1.75 * N0^0.75) * 1e18
//!   factor = 720 (sixth moment of exponential DSD)
//!   rho_x  = hydrometeor material density (kg m^-3)
//!   N0     = intercept parameter (m^-4)
//!   1e18   = unit conversion from m^3 to mm^6 m^-3
//!
//! Ice species (snow, graupel) are further scaled by |K_ice|^2/|K_water|^2
//! = 0.224 for the dielectric factor.

use rayon::prelude::*;

use crate::compute::ComputeOpts;
use crate::error::WrfResult;
use crate::file::WrfFile;

const RD: f64 = 287.058;

// --- Precomputed reflectivity coefficients (mm^6 m^-3) ---
//
// C_rain = 720 / (pi^1.75 * 1000^1.75 * (8e6)^0.75) * 1e18
// C_snow = 0.224 * 720 / (pi^1.75 * 100^1.75 * (2e7)^0.75) * 1e18
// C_graup = 0.224 * 720 / (pi^1.75 * 400^1.75 * (4e6)^0.75) * 1e18
//
// N0_rain  = 8.0e6  m^-4,  rho_water  = 1000 kg/m^3
// N0_snow  = 2.0e7  m^-4,  rho_snow   =  100 kg/m^3
// N0_graup = 4.0e6  m^-4,  rho_graupel=  400 kg/m^3
const COEFF_RAIN: f64 = 3.630_803_362_5e9;
const COEFF_SNOW: f64 = 2.300_359_648_2e10;
const COEFF_GRAUP: f64 = 6.798_580_734_2e9;

/// Simulated reflectivity (dBZ). `[nz, ny, nx]`
///
/// Matches WRF's `module_diag_functions.F` `REFL_10CM` diagnostic.
/// Rain, snow, and graupel contributions use Marshall-Palmer intercept
/// parameters and density corrections consistent with the WRF defaults
/// (WSM6/Thompson-like microphysics).
pub fn compute_dbz(f: &WrfFile, t: usize, _opts: &ComputeOpts) -> WrfResult<Vec<f64>> {
    let tk = f.temperature(t)?;
    let pres = f.full_pressure(t)?;
    let qv = f.qvapor(t)?;

    // Try to read hydrometeor fields; default to zero if absent
    let qr = f.read_var("QRAIN", t).unwrap_or_else(|_| vec![0.0; f.nxyz()]);
    let qs = f.read_var("QSNOW", t).unwrap_or_else(|_| vec![0.0; f.nxyz()]);
    let qg = f.read_var("QGRAUP", t).unwrap_or_else(|_| vec![0.0; f.nxyz()]);

    let n = f.nxyz();

    let dbz: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|i| {
            let t_k = tk[i];
            let qv_i = qv[i].max(0.0);
            let tv = t_k * (1.0 + 0.61 * qv_i);
            let rho = pres[i] / (RD * tv);

            let qr_i = qr[i].max(0.0);
            let qs_i = qs[i].max(0.0);
            let qg_i = qg[i].max(0.0);

            // Z_e (mm^6 m^-3) = C * (rho_air * q)^1.75
            let z_rain = if qr_i > 1e-8 {
                COEFF_RAIN * (rho * qr_i).powf(1.75)
            } else {
                0.0
            };

            let z_snow = if qs_i > 1e-8 {
                COEFF_SNOW * (rho * qs_i).powf(1.75)
            } else {
                0.0
            };

            let z_graup = if qg_i > 1e-8 {
                COEFF_GRAUP * (rho * qg_i).powf(1.75)
            } else {
                0.0
            };

            let z_total = z_rain + z_snow + z_graup;
            if z_total > 0.0 {
                10.0 * z_total.log10()
            } else {
                -30.0 // below noise floor
            }
        })
        .collect();

    Ok(dbz)
}

/// Maximum (composite) reflectivity (dBZ). `[ny, nx]`
pub fn compute_maxdbz(f: &WrfFile, t: usize, opts: &ComputeOpts) -> WrfResult<Vec<f64>> {
    let dbz_3d = compute_dbz(f, t, opts)?;
    let nxy = f.nxy();
    let nz = f.nz;

    let mut maxdbz = vec![f64::NEG_INFINITY; nxy];
    for k in 0..nz {
        let offset = k * nxy;
        for ij in 0..nxy {
            let val = dbz_3d[offset + ij];
            if val > maxdbz[ij] {
                maxdbz[ij] = val;
            }
        }
    }
    Ok(maxdbz)
}
