use check_ef_vasp::{read_from_dir, read_from_vaspout};
use clap::Parser;
use sigpipe;

use anyhow::Result;
use std::fmt::Write as _;
use std::fs::File;
use std::io::Write;

#[derive(Parser)]
#[command(
    about = "Check convergence of vasp geometry optimization",
    version,
    author
)]
pub struct Options {
    /// The directory where the calculation is performed.
    #[arg(short, long, default_value = ".")]
    pub dir: String,
    /// If turned on, forces are not checked.
    #[arg(short = 'f', long, action)]
    pub no_force: bool,
    /// If turned on, results are written to "convergence.dat".
    #[arg(short = 'w', long, action)]
    pub write_results: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    sigpipe::reset();
    let opts = Options::parse();

    // Try to read from vaspout.h5 and if it fails, read from OUTCAR
    let (energies, max_forces) = match read_from_vaspout(&opts.dir, !opts.no_force) {
        Ok((energies, max_forces)) => (energies, max_forces),
        Err(_) => read_from_dir(&opts.dir, !opts.no_force)?,
    };
    // Calculate relative energies
    let rel_e: &Vec<f64> = &energies.iter().map(|e| e - energies[0]).collect();
    let d_e: Vec<f64> = {
        let mut _padded_e = vec![0.0];
        _padded_e.extend(&energies);
        _padded_e[..].windows(2).map(|x| x[1] - x[0]).collect()
    };

    // Pretty print values
    let mut header = format!("{:<12}", "Step");
    match max_forces {
        Some(_) => write!(header, "{:<15}", "F_max (eV/A)")?,
        None => {}
    }
    write!(
        header,
        "{:<15}{:<15}{:<15}",
        "E_0 (eV)", "E - E0 (eV)", "dE (eV)"
    )?;

    println!("{}", header);
    println!("{}", vec!["-"; header.len()].join(""));
    let mut lines = vec![];
    for i in 0..energies.len() {
        let step = i + 1;
        let mut line = format!("{:<12}", step);
        if let Some(ref max_forces) = max_forces {
            write!(line, "{:<15.6}", max_forces[i])?;
        }
        write!(
            line,
            "{:<15.6}{:<15.6}{:<15.6}",
            energies[i], rel_e[i], d_e[i]
        )?;
        println!("{}", line);
        lines.push(line);
    }
    if opts.write_results {
        let mut file = File::create(std::path::Path::new(&opts.dir).join("convergence.dat"))?;
        file.write_all(header.as_bytes())?;
        file.write_all(b"\n")?;
        file.write_all(lines.join("\n").as_bytes())?;
    }
    Ok(())
}
