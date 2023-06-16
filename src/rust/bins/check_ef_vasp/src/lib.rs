use anyhow::Result;
use corelib::parser::outcar;
use std::fs::File;
use std::io::prelude::*;
use std::path::PathBuf;

pub fn read_file(filename: &PathBuf) -> Result<String> {
    let mut buf = String::new();
    let mut f = File::open(filename)?;
    f.read_to_string(&mut buf)?;
    Ok(buf)
}

pub fn read_max_forces(poscar: &str, outcar: &str) -> Vec<f64> {
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
