pub mod file_formats;
pub mod io;
pub mod parser;
pub mod periodic_table;
pub mod structure;

#[cfg(test)]
mod tests {
    use crate::structure::Structure;
    use nalgebra::{Const, Dyn, OMatrix, SMatrix};
    #[test]
    fn test_get_distance() {
        let test_structure = Structure {
            cell: 2.0 * SMatrix::<f64, 3, 3>::identity(),
            coords: OMatrix::<f64, Dyn, Const<3>>::from_row_slice(&[
                1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.9, 1.0, 1.0,
            ]),
            elements: vec!["H".to_string(), "O".to_string()],

            charge: 0.0,
        };
        println!("{:?}", test_structure.get_fractional_coords() * 2.0);
        println!("{:?}", test_structure.coords);
        println!("{:.3}", test_structure.get_distance(1, 2));
        assert_eq!(test_structure.get_distance(0, 1), (3.0_f64).sqrt());
    }
    #[test]
    fn test_periodic_table() {
        use crate::periodic_table;
        let result = periodic_table::get_atomic_number("H".to_string());
        assert_eq!(result, 1);
    }
}
