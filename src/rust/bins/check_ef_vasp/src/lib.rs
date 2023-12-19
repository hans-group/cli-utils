use anyhow;
use anyhow::Result;
use corelib::parser::oszicar;
use corelib::parser::outcar;
use ndarray_linalg::Norm;
use paris;
use std::fs::File;
use std::io::prelude::*;
use std::path::PathBuf;

use ndarray_stats::QuantileExt;

use hdf5;

use ndarray::prelude::*;
pub fn read_file(filename: &PathBuf) -> Result<String> {
    let mut buf = String::new();
    let mut f = File::open(filename)?;
    f.read_to_string(&mut buf)?;
    Ok(buf)
}

pub fn read_max_forces_outcar(poscar: &str, outcar: &str) -> Vec<f64> {
    outcar::read_forces(&poscar, &outcar)
        .iter()
        .map(|forces| {
            forces
                .iter()
                .map(|f| f[0].powi(2) + f[1].powi(2) + f[2].powi(2))
                .max_by(|&x, &y| x.total_cmp(&y))
                .unwrap()
                .sqrt()
        })
        .collect()
}

pub fn read_from_dir(
    dir: &str,
    read_force: bool,
) -> anyhow::Result<(Vec<f64>, Option<Vec<f64>>), anyhow::Error> {
    // Read all files to strings
    let mut logger = paris::Logger::new();
    let dir = std::path::Path::new(dir);
    let poscar = read_file(&dir.join("POSCAR"))?;
    let oszicar = read_file(&dir.join("OSZICAR"))?;
    let outcar = match read_force {
        true => Some(read_file(&dir.join("OUTCAR"))?),
        false => None,
    };

    // Read energies from oszicar
    let energies = oszicar::read_energies(&oszicar);
    if energies.is_empty() {
        logger.error("No SCF loop found.");
        std::process::exit(1);
    }
    // Calculate max forces
    let max_forces = match read_force {
        true => Some(read_max_forces_outcar(&poscar, &outcar.unwrap())),
        false => None,
    };

    Ok((energies, max_forces))
}

pub fn read_from_vaspout(
    dir: &str,
    read_force: bool,
) -> anyhow::Result<(Vec<f64>, Option<Vec<f64>>), anyhow::Error> {
    let filepath = std::path::Path::new(dir).join("vaspout.h5");
    let file = hdf5::File::open(filepath)?;
    let dynamics = file.group("intermediate/ion_dynamics")?;
    let energies: Vec<f64> = {
        let ds = dynamics.dataset("energies")?;
        let data = ds.read_2d::<f64>()?;
        let energies = data.column(2);
        // Convert to Vec<f64>
        energies.iter().map(|x| *x).collect()
    };

    // Force should be masked because force component of fixed dimension is considered
    // as zero in geometry optimization.
    // F_masked = F * mask (F, mask are both N x 3 matrix)
    let force_mask = {
        let poscar_group = file.group("input/poscar")?;
        let n_atoms = poscar_group
            .dataset("position_ions")?
            .read_2d::<f64>()?
            .shape()[0];

        let selective_dynamics_on = poscar_group
            .dataset("selective_dynamics")?
            .read_scalar::<i32>()?
            == 1;
        match selective_dynamics_on {
            true => poscar_group
                .dataset("selective_dynamics_ions")?
                .read_2d::<i32>()?
                .mapv(|x| x as f64),
            false => Array2::from_elem((n_atoms, 3), 1.0),
        }
    };
    let max_forces: Option<Vec<f64>> = match read_force {
        true => {
            let ds = dynamics.dataset("forces")?;
            let mut forces = ds.read_dyn::<f64>()?;
            let max_forces = forces
                .axis_iter_mut(Axis(0))
                .map(|x| {
                    let masked_force = x.to_owned() * &force_mask;
                    let norms = masked_force.map_axis(Axis(1), |x| x.norm_l2());
                    *norms.max().unwrap()
                })
                .collect();
            Some(max_forces)
        }
        false => None,
    };

    Ok((energies, max_forces))
}
