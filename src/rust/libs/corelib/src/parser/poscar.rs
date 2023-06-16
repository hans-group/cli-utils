pub fn n_atoms_in_poscar(poscar: &str) -> usize {
    poscar
        .lines()
        .nth(6)
        .unwrap()
        .split_whitespace()
        .map(|x| x.parse::<usize>().unwrap())
        .sum()
}

pub fn read_unit_cell(poscar: &str) -> Vec<Vec<f64>> {
    poscar
        .lines()
        .skip(2)
        .take(3)
        .map(|line| {
            line.trim()
                .split_whitespace()
                .map(|x| x.parse::<f64>().unwrap())
                .collect()
        })
        .collect()
}

pub fn get_fix_mask(poscar: &str, n_atoms: usize) -> Vec<f64> {
    match poscar.lines().nth(7).unwrap().trim().starts_with('S') {
        true => poscar
            .lines()
            .skip(9)
            .take(n_atoms)
            .map(|line| match line.contains("F") {
                true => 0.0,
                false => 1.0,
            })
            .collect(),
        false => vec![1.0; n_atoms],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const POSCAR_STR: &str = "H2O molecule
    1.0
    13.00  0.00  0.00
    0.00  14.00  0.00
    0.00  0.00  15.00
    O   H
    1   2
    Selective dynamics
    Cartesian
    6.500000  7.000000  7.798154 F F F
    6.500000  7.763240  7.201845 T T T
    6.500000  6.236761  7.201845 T T T
";

    #[test]
    fn test_n_atoms_in_poscar() {
        let n_atoms = n_atoms_in_poscar(POSCAR_STR);
        assert_eq!(n_atoms, 3);
    }

    #[test]
    fn test_read_unit_cell() {
        let unit_cell = read_unit_cell(POSCAR_STR);
        assert_eq!(
            unit_cell,
            vec![
                vec![13.0, 0.0, 0.0],
                vec![0.0, 14.0, 0.0],
                vec![0.0, 0.0, 15.0]
            ]
        );
    }

    #[test]
    fn test_get_fix_mask() {
        let fix_mask = get_fix_mask(POSCAR_STR, 3);
        assert_eq!(fix_mask, vec![0.0, 1.0, 1.0]);
    }
}
