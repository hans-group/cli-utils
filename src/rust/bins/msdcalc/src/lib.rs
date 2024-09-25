extern crate ndarray;
extern crate ndarray_linalg;

use ndarray::prelude::*;
use ndarray_linalg::Inverse;

pub fn apply_mic(dr: &Array1<f64>, cell: &Array2<f64>, inv_cell: &Array2<f64>) -> Array1<f64> {
    let dr_mic = {
        let s = dr.dot(inv_cell);
        let s = s.mapv(|x| x - x.round());
        s.dot(cell)
    };
    dr_mic
}

pub fn unwrap_positions(
    positions: &Array3<f64>,
    cell: &Array2<f64>,
    inv_cell: &Array2<f64>,
) -> Array3<f64> {
    let (n_frames, n_atoms, _) = positions.dim();
    let mut unwrapped: Array3<f64> = Array3::zeros((n_frames, n_atoms, 3));
    unwrapped
        .slice_mut(s![0, .., ..])
        .assign(&positions.slice(s![0, .., ..]));
    for i in 1..n_frames {
        let mut diff = &positions.slice(s![i, .., ..]) - &positions.slice(s![i - 1, .., ..]);
        for j in 0..n_atoms {
            let diff_j = apply_mic(&diff.slice(s![j, ..]).to_owned(), cell, inv_cell);
            diff.slice_mut(s![j, ..]).assign(&diff_j);
        }
        let update = &unwrapped.slice(s![i - 1, .., ..]) + &diff;
        unwrapped.slice_mut(s![i, .., ..]).assign(&update);
    }

    unwrapped
}

fn calculate_msd_array(unwrapped: &Array3<f64>, max_delta: usize) -> Array1<f64> {
    let (n_frames, n_atoms, _) = unwrapped.dim();
    let mut msd: Array1<f64> = Array1::zeros(max_delta);
    let mut count: Array1<f64> = Array1::zeros(max_delta);
    for i in 0..n_frames {
        for j in i + 1..(i + max_delta + 1).min(n_frames) - 1 {
            let dt = j - i;
            let disp = &unwrapped.slice(s![j, .., ..]) - &unwrapped.slice(s![i, .., ..]);
            msd[dt] += &disp.mapv(|x| x * x).sum();
            count[dt] += n_atoms as f64;
        }
    }

    msd / &count
}

pub fn calculate_msd(
    positions: Array3<f64>,
    cell: &Array2<f64>,
    elements: Vec<String>,
    selected_elements: Option<Vec<String>>,
    max_time_delta: Option<f64>,
    dt: f64,
    stride: usize,
) -> (Array1<f64>, Array1<f64>) {
    let n_frames = positions.dim().0;
    let n_atoms = positions.dim().1;
    let inv_cell = cell.inv().unwrap();

    let atom_indices = match selected_elements {
        Some(selected_elements) => {
            let mut atom_indices = vec![];
            for element in selected_elements {
                let indices = elements
                    .iter()
                    .enumerate()
                    .filter(|(_, e)| **e == element)
                    .map(|(i, _)| i)
                    .collect::<Vec<usize>>();
                atom_indices.extend(indices);
            }
            atom_indices
        }
        None => (0..n_atoms).collect(),
    };
    let atom_indices = Array1::<usize>::from(atom_indices);
    let unwrapped_positions = unwrap_positions(&positions, cell, &inv_cell);
    let max_delta = match max_time_delta {
        Some(max_time_delta) => ((max_time_delta / (dt * stride as f64)) as usize).min(n_frames - 1),
        None => n_frames - 1,
    };
    println!("Calculating MSD for {} frames: total {:.1} ps", max_delta, max_time_delta.unwrap_or(n_frames as f64 * dt * stride as f64));
    let elem_positions = {
        let mut pos = Array3::zeros((n_frames, atom_indices.len(), 3));
        for (ii, &i) in atom_indices.iter().enumerate() {
            pos.slice_mut(s![.., ii, ..])
                .assign(&unwrapped_positions.slice(s![.., i, ..]));
        }
        pos
    };
    let mut msd = calculate_msd_array(&elem_positions, max_delta);
    msd[0] = 0.0;
    let times = Array1::range(0.0, max_delta as f64, 1.0) * dt * stride as f64;

    (times, msd)
}
