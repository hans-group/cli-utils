use super::poscar::{get_fix_mask, n_atoms_in_poscar};
use rayon::prelude::*;

enum BlockProperty {
    Forces,
    Positions,
}

fn read_frame_block(block: &str, mask: &[f64], property: BlockProperty) -> Vec<Vec<f64>> {
    let skip_size = match property {
        BlockProperty::Positions => 0,
        BlockProperty::Forces => 3,
    };
    block
        .lines()
        .enumerate()
        .map(|(i, line)| {
            line.split_ascii_whitespace()
                .skip(skip_size)
                .take(3)
                .map(|x| mask[i] * x.parse::<f64>().unwrap())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>()
}

pub fn read_forces(poscar: &str, outcar: &str) -> Vec<Vec<Vec<f64>>> {
    let n_atoms: usize = n_atoms_in_poscar(poscar);
    let mask = get_fix_mask(&poscar, n_atoms);
    let outcar_lines = outcar.lines().collect::<Vec<_>>();

    outcar_lines
        .par_iter()
        .enumerate()
        .filter(|(_, line)| line.contains("TOTAL-FORCE"))
        .map(|(i, _)| i)
        .map(|i| {
            let block = &outcar_lines[i + 2..i + 2 + n_atoms].join("\n");
            read_frame_block(&block, &mask, BlockProperty::Forces)
        })
        .collect()
}

pub fn read_trajectory(poscar: &str, outcar: &str) -> Vec<Vec<Vec<f64>>> {
    let n_atoms: usize = n_atoms_in_poscar(poscar);
    let outcar_lines = outcar.lines().collect::<Vec<_>>();
    let ones = vec![1.0; n_atoms];

    outcar_lines
        .par_iter()
        .enumerate()
        .filter(|(_, line)| line.contains("TOTAL-FORCE"))
        .map(|(i, _)| i)
        .map(|i| {
            let block = &outcar_lines[i + 2..i + 2 + n_atoms].join("\n");
            read_frame_block(&block, &ones, BlockProperty::Positions)
        })
        .collect()
}
