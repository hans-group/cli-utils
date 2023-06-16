use nalgebra::Vector3;
use nalgebra::{Matrix3, RowVector3};
use std::fmt::Display;
use std::fs::read_to_string;
use std::str::FromStr;
use thiserror::Error;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Lattice {
    pub a: [f64; 3],
    pub b: [f64; 3],
    pub c: [f64; 3],
}

impl Lattice {
    pub fn new(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Self {
        Lattice { a, b, c }
    }
}

impl Lattice {
    pub fn lattice_params(self) -> LatticeParams {
        self.into()
    }

    pub fn to_matrix(self) -> Matrix3<f64> {
        let a_vec = RowVector3::new(self.a[0], self.a[1], self.a[2]);
        let b_vec = RowVector3::new(self.b[0], self.b[1], self.b[2]);
        let c_vec = RowVector3::new(self.c[0], self.c[1], self.c[2]);
        Matrix3::from_rows(&[a_vec, b_vec, c_vec])
    }

    pub fn to_vec(self) -> Vec<Vec<f64>> {
        vec![self.a.to_vec(), self.b.to_vec(), self.c.to_vec()]
    }

    pub fn volume(self) -> f64 {
        self.to_matrix().determinant().abs()
    }

    pub fn reciprocal(self) -> Lattice {
        let pi = std::f64::consts::PI;
        let volume = self.volume();

        let mat = self.to_matrix();
        let a_vec = &mat.row(0);
        let b_vec = &mat.row(1);
        let c_vec = &mat.row(2);

        let a_reciprocal = 2.0 * pi / volume * b_vec.cross(&c_vec);
        let b_reciprocal = 2.0 * pi / volume * c_vec.cross(&a_vec);
        let c_reciprocal = 2.0 * pi / volume * a_vec.cross(&b_vec);
        Lattice {
            a: [a_reciprocal[0], a_reciprocal[1], a_reciprocal[2]],
            b: [b_reciprocal[0], b_reciprocal[1], b_reciprocal[2]],
            c: [c_reciprocal[0], c_reciprocal[1], c_reciprocal[2]],
        }
    }
}

impl Display for Lattice {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{:.9}  {:.9}  {:.9}\n{:.9}  {:.9}  {:.9}\n{:.9}  {:.9}  {:.9}",
            self.a[0],
            self.a[1],
            self.a[2],
            self.b[0],
            self.b[1],
            self.b[2],
            self.c[0],
            self.c[1],
            self.c[2]
        )
    }
}

impl From<Lattice> for LatticeParams {
    fn from(lattice: Lattice) -> Self {
        let a_vec = Vector3::new(lattice.a[0], lattice.a[1], lattice.a[2]);
        let b_vec = Vector3::new(lattice.b[0], lattice.b[1], lattice.b[2]);
        let c_vec = Vector3::new(lattice.c[0], lattice.c[1], lattice.c[2]);
        let a = a_vec.norm();
        let b = b_vec.norm();
        let c = c_vec.norm();

        let alpha = b_vec.angle(&c_vec).to_degrees();
        let beta = a_vec.angle(&c_vec).to_degrees();
        let gamma = a_vec.angle(&b_vec).to_degrees();
        LatticeParams::new(a, b, c, alpha, beta, gamma)
    }
}

#[derive(Debug)]
pub struct LatticeParams {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub alpha: f64,
    pub beta: f64,
    pub gamma: f64,
}

impl LatticeParams {
    #[rustfmt::skip]
    pub fn new(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> Self {
        LatticeParams { a, b, c, alpha, beta, gamma }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordinateSystem {
    Cartesian,
    Direct,
}

impl FromStr for CoordinateSystem {
    type Err = PoscarParseError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            s if s.to_ascii_lowercase().starts_with('c') => Ok(CoordinateSystem::Cartesian),
            s if s.to_ascii_lowercase().starts_with('d') => Ok(CoordinateSystem::Direct),
            _ => Err(PoscarParseError::UnknownCoordinateSystem(s.to_string())),
        }
    }
}

impl Display for CoordinateSystem {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            CoordinateSystem::Cartesian => write!(f, "Cartesian"),
            CoordinateSystem::Direct => write!(f, "Direct"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SelectiveDynamics(bool);

impl FromStr for SelectiveDynamics {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "T" => Ok(SelectiveDynamics(true)),
            "F" => Ok(SelectiveDynamics(false)),
            _ => Err(format!("Invalid selective dynamics: {}", s)),
        }
    }
}

impl Display for SelectiveDynamics {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", if self.0 { "T" } else { "F" })
    }
}

#[derive(Debug)]
pub struct Poscar {
    pub lattice: Lattice,
    pub species: Vec<String>,
    pub num_species: Vec<usize>,
    pub chemical_symbols: Vec<String>,
    pub positions: Vec<[f64; 3]>,
    pub coord_type: CoordinateSystem,
    pub selective_dynamics: Option<Vec<[SelectiveDynamics; 3]>>,
}

#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum PoscarParseError {
    #[error("File not found: {0}")]
    FileNotFound(String),
    #[error("Unknown coordinate system {0}.")]
    UnknownCoordinateSystem(String),
    #[error("Unknown selective dynamics {0}.")]
    BadSelectiveDynamics(String),
    #[error("Malformed POSCAR file.")]
    BadPOSCAR,
    #[error("BadScalingFactor")]
    BadScalingFactor,
    #[error("BadLattice")]
    BadLattice,
}

use PoscarParseError::*;
impl FromStr for Poscar {
    type Err = PoscarParseError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut lines = s.lines();
        lines.next().ok_or(BadPOSCAR)?;

        let scaling_factor = lines
            .next()
            .ok_or(BadScalingFactor)?
            .trim()
            .parse::<f64>()
            .map_err(|_| BadScalingFactor)?;

        let lattice = {
            let mut lattice = [[0.0; 3]; 3];
            for i in 0..3 {
                let line = lines.next().ok_or(BadLattice)?;
                let mut tokens = line.trim().split_whitespace();
                for j in 0..3 {
                    lattice[i][j] = tokens
                        .next()
                        .ok_or(BadLattice)?
                        .parse::<f64>()
                        .map_err(|_| BadLattice)?
                        * scaling_factor;
                }
            }
            Lattice::new(lattice[0], lattice[1], lattice[2])
        };

        let species = lines
            .next()
            .ok_or(BadPOSCAR)?
            .trim()
            .split_whitespace()
            .map(|s| s.to_string())
            .collect::<Vec<String>>();
        let num_species = lines
            .next()
            .ok_or(BadPOSCAR)?
            .trim()
            .split_whitespace()
            .map(|s| s.parse::<usize>().map_err(|_| BadPOSCAR))
            .collect::<Result<Vec<usize>, PoscarParseError>>()?;
        let chemical_symbols = species
            .iter()
            .zip(&num_species)
            .flat_map(|(s, n)| std::iter::repeat(s.clone()).take(*n))
            .collect::<Vec<_>>();

        let coord_or_selective_dynamics = lines.next().ok_or(BadPOSCAR)?;
        let mut coord_type = coord_or_selective_dynamics.to_string().trim().to_string();
        let is_selective_dynamics = coord_or_selective_dynamics.trim().starts_with("S");
        if is_selective_dynamics {
            coord_type = lines.next().ok_or(BadPOSCAR)?.trim().to_string();
        }

        let mut positions = Vec::new();
        let mut selective_dynamics = Vec::new();
        for _ in 0..num_species.iter().sum::<usize>() {
            let line = lines.next().ok_or(BadPOSCAR)?.trim();
            let mut tokens = line.split_whitespace();
            let x = tokens
                .next()
                .ok_or(BadPOSCAR)?
                .trim()
                .parse::<f64>()
                .map_err(|_| BadPOSCAR)?;
            let y = tokens
                .next()
                .ok_or(BadPOSCAR)?
                .trim()
                .parse::<f64>()
                .map_err(|_| BadPOSCAR)?;
            let z = tokens
                .next()
                .ok_or(BadPOSCAR)?
                .trim()
                .parse::<f64>()
                .map_err(|_| BadPOSCAR)?;
            positions.push([x, y, z]);
            if is_selective_dynamics {
                let fix_x: SelectiveDynamics = tokens
                    .next()
                    .ok_or(BadPOSCAR)?
                    .trim()
                    .parse()
                    .map_err(|s| BadSelectiveDynamics(s))?;
                let fix_y: SelectiveDynamics = tokens
                    .next()
                    .ok_or(BadPOSCAR)?
                    .trim()
                    .parse()
                    .map_err(|s| BadSelectiveDynamics(s))?;
                let fix_z: SelectiveDynamics = tokens
                    .next()
                    .ok_or(BadPOSCAR)?
                    .trim()
                    .parse()
                    .map_err(|s| BadSelectiveDynamics(s))?;
                selective_dynamics.push([fix_x, fix_y, fix_z]);
            }
        }

        Ok(Poscar {
            lattice,
            species,
            num_species,
            chemical_symbols,
            positions,
            coord_type: coord_type.parse()?,
            selective_dynamics: if !selective_dynamics.is_empty() {
                Some(selective_dynamics)
            } else {
                None
            },
        })
    }
}

impl Display for Poscar {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "Poscar")?;
        writeln!(f, "1.0")?;
        writeln!(f, "{}", self.lattice)?;
        writeln!(f, "{}", self.species.join(" "))?;
        writeln!(
            f,
            "{}",
            self.num_species
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(" ")
        )?;
        if self.selective_dynamics.is_some() {
            writeln!(f, "Selective dynamics")?;
        }
        writeln!(f, "{}", self.coord_type)?;
        for i in 0..self.positions.len() {
            write!(
                f,
                "{:.9}  {:.9}  {:.9}  ",
                self.positions[i][0], self.positions[i][1], self.positions[i][2]
            )?;
            if let Some(sd) = &self.selective_dynamics {
                write!(f, "{}  {}  {}  ", sd[i][0], sd[i][1], sd[i][2])?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

impl Poscar {
    pub fn from_file(path: &str) -> Result<Self, PoscarParseError> {
        let contents = read_to_string(path).map_err(|_| FileNotFound(path.to_string()))?;
        Poscar::from_str(&contents)
    }

    pub fn get_lattice_params(&self) -> LatticeParams {
        self.lattice.into()
    }

    pub fn get_reciprocal_lattice(&self) -> Lattice {
        self.lattice.reciprocal()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const POSCAR: &str = "lkjh
    1.0000000000000000
        5.5437171645025325    0.0000000000000000    0.0000000000000000
        2.7718585822512662    4.8009998958550284    0.0000000000000000
        0.0000000000000000    0.0000000000000000   24.0528522208933317
    Pt
    3
    Selective dynamics
    Direct
    -0.7384622484064965 -1.2288412475232857  0.0000000000000000   F   F   F
    -0.2384622484064966 -1.2288412475232857  0.0000000000000000   F   F   F
    -0.7384622484064967 -0.7288412475232857  0.0000000000000000   F   F   F
";

    #[test]
    fn test_lattice_params() {
        let lattice = Lattice {
            a: [1.0, 0.0, 0.0],
            b: [0.0, 1.0, 0.0],
            c: [0.0, 0.0, 1.0],
        };
        let lattice_params = LatticeParams::from(lattice);
        assert_eq!(lattice_params.a, 1.0);
        assert_eq!(lattice_params.b, 1.0);
        assert_eq!(lattice_params.c, 1.0);
        assert_eq!(lattice_params.alpha, 90.0);
        assert_eq!(lattice_params.beta, 90.0);
        assert_eq!(lattice_params.gamma, 90.0);
    }

    #[test]
    fn test_reciprocal_lattice() {
        let pi = std::f64::consts::PI;
        let poscar = Poscar::from_str(POSCAR).unwrap();
        let reciprocal_lattice = poscar.get_reciprocal_lattice();

        assert!((reciprocal_lattice.a[0] / (2.0 * pi) - 1.80384383e-01) < 1e-6);
        assert!((reciprocal_lattice.a[1] / (2.0 * pi) - -1.04144972e-01) < 1e-6);
        assert!((reciprocal_lattice.a[2] / (2.0 * pi) - 0.0) < 1e-6);

        assert!((reciprocal_lattice.b[0] / (2.0 * pi) - 0.0) < 1e-6);
        assert!((reciprocal_lattice.b[1] / (2.0 * pi) - 2.08289944e-01) < 1e-6);
        assert!((reciprocal_lattice.b[2] / (2.0 * pi) - 0.0) < 1e-6);

        assert!((reciprocal_lattice.c[0] / (2.0 * pi) - 0.0) < 1e-6);
        assert!((reciprocal_lattice.c[1] / (2.0 * pi) - 0.0) < 1e-6);
        assert!((reciprocal_lattice.c[2] / (2.0 * pi) - 4.15751110e-02) < 1e-6);

        println!("{:?}", reciprocal_lattice.lattice_params());
    }
}
