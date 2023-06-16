use std::fmt::Display;

#[derive(Debug)]
pub struct XYZBlock {
    pub elements: Vec<String>,
    pub positions: Vec<Vec<f64>>,
    pub cell: Option<Vec<Vec<f64>>>,
}

fn format_cell_as_str(cell: &Vec<Vec<f64>>) -> String {
    format!(
        "{:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6}",
        cell[0][0],
        cell[0][1],
        cell[0][2],
        cell[1][0],
        cell[1][1],
        cell[1][2],
        cell[2][0],
        cell[2][1],
        cell[2][2]
    )
}

impl Display for XYZBlock {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let prop_header = "Properties=species:S:1:pos:R:3";
        // Atom string
        let num_atoms = self.elements.len();
        let xyz_header = match &self.cell {
            Some(cell) => {
                let cell_str = format_cell_as_str(cell);
                format!(
                    "{}\nLattice=\"{}\" {} pbc=\"T T T\"",
                    num_atoms, cell_str, prop_header
                )
            }
            None => format!("{}\n{}", num_atoms, prop_header),
        };

        write!(f, "{}", xyz_header)?;
        for (i, element) in self.elements.iter().enumerate() {
            let x = self.positions[i][0];
            let y = self.positions[i][1];
            let z = self.positions[i][2];
            write!(
                f,
                "\n{:<3}   {:<15.9}  {:<15.9}  {:<15.9}",
                element, x, y, z
            )?;
        }
        Ok(())
    }
}
