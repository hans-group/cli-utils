use std::io::Write;

use clap::Parser;
use corelib::file_formats::poscar::{Lattice, Poscar};
use spglib::cell::Cell;

#[derive(Parser)]
#[command(
    version = "0.1.0",
    author = "Minjoon Hong",
    about = "Symmetrize POSCAR file"
)]
struct Opts {
    /// Input POSCAR file
    #[arg(short, long)]
    input: String,
    /// Tolernace for symmetry detection
    #[arg(short, long, default_value = "1e-6")]
    symprec: f64,
}

fn main() {
    let args = Opts::parse();
    let filename = args.input;
    let mut poscar = Poscar::from_file(&filename).unwrap();

    let lattice = [
        poscar.lattice.a.clone(),
        poscar.lattice.b.clone(),
        poscar.lattice.c.clone(),
    ];
    let positions = poscar.positions.clone();
    let types = poscar.chemical_symbols.clone();

    let mut typemap = std::collections::HashMap::new();
    let mut typenums = vec![];
    for (i, t) in poscar.species.iter().enumerate() {
        typemap.insert(t, i);
    }
    for t in types {
        typenums.push(typemap[&t] as i32);
    }
    println!("typemap: {:?}", typemap);
    println!("typenums: {:?}", typenums);

    let mut cell = Cell::new(&lattice, &positions, &typenums);
    cell.standardize(true, false, args.symprec).unwrap();
    poscar.positions = cell.positions;
    poscar.lattice = Lattice::new(cell.lattice[0], cell.lattice[1], cell.lattice[2]);

    let poscar_str = poscar.to_string();
    //write to file
    let mut file = std::fs::File::create("symmetrized.vasp").unwrap();
    file.write_all(poscar_str.as_bytes()).unwrap();
    println!("\nDone.symmetrized.vasp written");
}
