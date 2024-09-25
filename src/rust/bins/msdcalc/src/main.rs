use anyhow;
use clap::Parser;
use extxyz::{read_xyz_frames, Info, RawAtoms};
use msdcalc::calculate_msd;
use ndarray::prelude::*;
use std::io::Write;

#[derive(Parser)]
#[command(about = "Calculate the mean square displacement", version, author)]
pub struct Options {
    /// The input (ext)xyz file that contains the MD trajectory.
    #[arg(short, long)]
    pub input: String,
    /// The output file where the MSD is written.
    #[arg(short, long)]
    pub output: String,
    /// elements to calculate MSD for.
    #[arg(short, long, value_delimiter = ',')]
    pub species: Vec<String>,
    /// The number of frames to skip.
    #[arg(long, default_value = "0")]
    pub skip_frames: usize,
    /// The interval of frames
    #[arg(short = 'S', long, default_value = "1")]
    pub stride: usize,
    /// The time interval between frames in ps.
    #[arg(long, default_value = "0.001")]
    pub dt: f64,
    /// Maximum time difference to consider. If None, all frames are considered.
    #[arg(short = 'M', long)]
    pub max_time_delta: Option<f64>,
}

fn main() {
    let opts = Options::parse();
    let frames = read_extxyz(&opts.input, Some(opts.skip_frames), None, Some(opts.stride)).unwrap();
    let cell = match &frames[0].cell {
        Some(cell) => cell.clone(),
        None => panic!("No cell information found in the input file"),
    };
    let positions = {
        let mut positions: Array3<f64> = Array3::zeros((frames.len(), frames[0].symbols.len(), 3));
        for (i, frame) in frames.iter().enumerate() {
            positions.slice_mut(s![i, .., ..]).assign(&frame.positions);
        }
        positions
    };
    let symbols = frames[0].symbols.clone();
    println!("Number of selected frames: {}", positions.dim().0);
    println!("Calculating MSD for species: {:?}", opts.species);
    let (times, msd) = calculate_msd(
        positions,
        &cell,
        symbols,
        Some(opts.species),
        opts.max_time_delta,
        opts.dt,
        opts.stride,
    );
    let mut file = std::fs::File::create(&opts.output).unwrap();
    writeln!(file, "Time(ps) MSD (A^2)").unwrap();
    for (x, y) in times.iter().zip(msd.iter()) {
        writeln!(file, "{:.6} {:.6}", x, y).unwrap();
    }
}

pub struct Atoms {
    pub symbols: Vec<String>,
    pub positions: Array2<f64>,
    pub cell: Option<Array2<f64>>,
}

pub fn read_extxyz(
    filename: impl AsRef<std::path::Path> + std::fmt::Display,
    start: Option<usize>,
    end: Option<usize>,
    step: Option<usize>,
) -> anyhow::Result<Vec<Atoms>> {
    let selection = match (start, end, step) {
        (Some(start), Some(end), Some(step)) => (start..end).step_by(step),
        (Some(start), Some(end), None) => (start..end).step_by(1),
        (Some(start), None, None) => (start..usize::MAX).step_by(1),
        (Some(start), None, Some(step)) => (start..usize::MAX).step_by(step),
        (None, Some(end), Some(step)) => (0..end).step_by(step),
        (None, Some(end), None) => (0..end).step_by(1),
        (None, None, None) => (0..usize::MAX).step_by(1),
        (None, None, Some(step)) => (0..usize::MAX).step_by(step),
    };
    let frames = read_xyz_frames(filename, selection)?;

    let mut atoms_list = vec![];
    for frame in frames {
        let atoms = RawAtoms::parse_from(&frame)?;
        // it will returen error if the comment is not in normal extxyz format
        let info: Info = atoms.comment.parse()?;
        // get molecule's properties
        let lattice = {
            let raw = info.get("Lattice").map(|x| x.as_array().unwrap());
            match raw {
                Some(vals) => {
                    let vals = vals.as_slice();
                    let mut lattice = Array2::zeros((3, 3));
                    for i in 0..3 {
                        for j in 0..3 {
                            lattice[[i, j]] = vals[i * 3 + j].as_f64().unwrap();
                        }
                    }
                    Some(lattice)
                }
                None => None,
            }
        };

        // get atom's properties
        let mut symbols = vec![];
        let mut positions = vec![];

        for atom in atoms.atoms {
            symbols.push(atom.element.to_string());
            positions.extend_from_slice(atom.position.as_slice());
        }
        let positions = Array2::from_shape_vec((symbols.len(), 3), positions)?;

        // Parse force

        let atoms = Atoms {
            symbols,
            positions,
            cell: lattice,
        };
        atoms_list.push(atoms);
    }

    Ok(atoms_list)
}
