use nalgebra::{Const, Dyn, OMatrix, SMatrix};

use crate::periodic_table::ATOMIC_NUMBERS;
pub struct Structure {
    pub cell: SMatrix<f64, 3, 3>,
    // matrix with known number of columns, but unknown number of rows
    pub coords: OMatrix<f64, Dyn, Const<3>>,
    pub elements: Vec<String>,
    pub charge: f64,
}

impl Structure {
    pub fn get_fractional_coords(&self) -> OMatrix<f64, Dyn, Const<3>> {
        let inv_cell = self.cell.try_inverse().unwrap();
        (inv_cell * self.coords.transpose()).transpose()
    }

    pub fn get_atomic_numbers(&self) -> Vec<u32> {
        self.elements.iter().map(|e| ATOMIC_NUMBERS[e]).collect()
    }

    pub fn get_distance(&self, i: usize, j: usize) -> f64 {
        let displacement = self.periodic_displacement(i, j);
        displacement.norm()
    }

    // angle is defined as ijk
    pub fn get_angle(&self, i: usize, j: usize, k: usize) -> f64 {
        let displacement_ij = self.periodic_displacement(i, j);
        let displacement_jk = self.periodic_displacement(j, k);
        displacement_ij.angle(&displacement_jk)
    }

    fn periodic_displacement(&self, i: usize, j: usize) -> OMatrix<f64, Const<1>, Const<3>> {
        let frac_coords = self.get_fractional_coords();
        let displacement = frac_coords.row(j) - frac_coords.row(i);
        let periodic_disp = displacement.map(|x| (x + 0.5) % 1.0 - 0.5);

        periodic_disp * self.cell
    }
}
