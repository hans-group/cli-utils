use clap::Parser;
use corelib::file_formats::kpoints::{Kpoints, KpointsGridScheme};
use corelib::file_formats::poscar::Poscar;
use std::fs::File;
use std::io::Write;

const TARGET_FILLEPATH: &str = "KPOINTS";

/// Generates KPOINTS file for VASP calculations.
#[derive(Parser, Debug)]
#[command(author,version,about,long_about = None)]
struct Args {
    /// The scheme to generate the k-points grid.
    /// Defaults to Gamma, which means the program will generate a Gamma-centered grid.
    /// Available options are:
    /// - Gamma (or G, g): Gamma-centered grid
    /// - Monkhorst (or M, m): Monkhorst-Pack grid
    #[arg(short, long, value_parser=parse_grid_scheme, default_value = "Gamma")]
    grid_scheme: KpointsGridScheme,
    /// The number of k-points along each direction.
    /// Defaults to 1 1 1, which means the program will generate a single k-point.
    #[arg(short, long, value_delimiter = ' ', num_args = 3)]
    ngrids: Option<Vec<u32>>,
    /// The density of k-points along each direction.
    /// Requires POSCAR file.
    #[arg(short, long)]
    density: Option<f64>,
}

fn parse_grid_scheme(s: &str) -> Result<KpointsGridScheme, String> {
    match s {
        s if s.to_ascii_lowercase().starts_with('g') => Ok(KpointsGridScheme::Gamma),
        s if s.to_ascii_lowercase().starts_with('m') => Ok(KpointsGridScheme::MonkhorstPack),
        _ => Err(format!("Unknown scheme {}.", s)),
    }
}

fn main() {
    let args = Args::parse();
    // Write the k-points grid to the file.
    match (args.ngrids, args.density) {
        (Some(_), Some(_)) => {
            println!("Cannot specify both ngrids and density.");
            std::process::exit(1);
        }
        (None, None) => {
            println!("Must specify either ngrids or density.");
            std::process::exit(1);
        }
        (Some(ngrids), None) => {
            let kpoints = Kpoints::new(args.grid_scheme, [ngrids[0], ngrids[1], ngrids[2]]);
            let mut file = File::create(TARGET_FILLEPATH).expect("Unable to create file");
            file.write_all(kpoints.to_string().as_bytes())
                .expect("Unable to write data");
        }
        (None, Some(density)) => {
            let poscar = Poscar::from_file("POSCAR").expect("Unable to read POSCAR file");
            let mut file = File::create(TARGET_FILLEPATH).expect("Unable to create file");
            let kpoints = Kpoints::from_density(args.grid_scheme, density, &poscar.lattice);
            file.write_all(kpoints.to_string().as_bytes())
                .expect("Unable to write data");
        }
    }
}
