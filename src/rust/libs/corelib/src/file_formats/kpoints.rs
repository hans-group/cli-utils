use crate::file_formats::poscar::Lattice;
use std::{fmt::Display, str::FromStr};
use thiserror::Error;

#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum KpointsError {
    #[error("Unknown scheme {0}.")]
    UnknownScheme(String),
    #[error("Zero mesh is not allow for k-points.")]
    ZeroMesh(String),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Kpoints {
    pub scheme: KpointsGridScheme,
    pub mesh: [u32; 3],
}

impl Display for Kpoints {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Autometic mesh\n 0\n{}\n{}  {}  {}\n0  0  0\n",
            self.scheme, self.mesh[0], self.mesh[1], self.mesh[2]
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KpointsGridScheme {
    Gamma,
    MonkhorstPack,
}

impl Display for KpointsGridScheme {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            KpointsGridScheme::Gamma => write!(f, "Gamma"),
            KpointsGridScheme::MonkhorstPack => write!(f, "Monkhorst-Pack"),
        }
    }
}

impl FromStr for KpointsGridScheme {
    type Err = KpointsError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            s if s.to_ascii_lowercase().starts_with('g') => Ok(KpointsGridScheme::Gamma),
            s if s.to_ascii_lowercase().starts_with('m') => Ok(KpointsGridScheme::MonkhorstPack),
            _ => Err(KpointsError::UnknownScheme(s.to_string())),
        }
    }
}

impl Kpoints {
    pub fn default() -> Kpoints {
        Kpoints {
            scheme: KpointsGridScheme::Gamma,
            mesh: [1, 1, 1],
        }
    }

    pub fn new(scheme: KpointsGridScheme, mesh: [u32; 3]) -> Result<Kpoints, KpointsError> {
        for m in mesh {
            if m == 0 {
                return Err(KpointsError::ZeroMesh(
                    "Zero mesh is not allow for k-points.".to_string(),
                ));
            }
        }
        Ok(Kpoints { scheme, mesh })
    }

    pub fn from_density(scheme: KpointsGridScheme, density: f64, lattice: &Lattice) -> Kpoints {
        let reciprocal_lattice_params = lattice.reciprocal().lattice_params();
        let a = reciprocal_lattice_params.a;
        let b = reciprocal_lattice_params.b;
        let c = reciprocal_lattice_params.c;

        let num_ka = std::cmp::max(1, (a * density).round() as u32); // Ensure at least 1 k-point
        let num_kb = std::cmp::max(1, (b * density).round() as u32);
        let num_kc = std::cmp::max(1, (c * density).round() as u32);
        let mesh = [num_ka, num_kb, num_kc];
        Kpoints { scheme, mesh }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::file_formats::poscar::Poscar;

    const POSCAR_TEST_BULK: &str = "Cu Fe Pt 
 1.0000000000000000
    30.0000000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000   30.0000000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000   30.0000000000000000
Cu  
3
Direct
0.3040000000000000  0.3040000000000000  0.5000000000000000
0.3040000000000000  0.3693333333333333  0.4346666666666666
0.3693333333333333  0.3040000000000000  0.4346666666666666
";
    #[test]
    fn test_kpoints() {
        let kpoints = Kpoints::new(KpointsGridScheme::Gamma, [3, 11, 7]).unwrap();
        println!("{}", kpoints);
        assert_eq!(
            kpoints.to_string(),
            "Autometic mesh\n 0\nGamma\n3  11  7\n0  0  0\n"
        );
    }

    #[test]
    fn test_kpoints_from_density() {
        let poscar = Poscar::from_str(POSCAR_TEST_BULK).unwrap();
        let kpoints = Kpoints::from_density(KpointsGridScheme::Gamma, 4.0, &poscar.lattice);
        println!("{}", kpoints);
    }
}
