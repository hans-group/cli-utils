use comfy_table::presets::UTF8_BORDERS_ONLY;
use comfy_table::*;
use std::process::Command;

const NODE_OWNERS: &str = include_str!("../data/hpc_node_owners.csv");

fn main() {
    let mut table = Table::new();

    let node_owners: Vec<&str> = NODE_OWNERS
        .lines()
        .skip(1)
        .map(|s| s.split(',').last().unwrap())
        .collect();

    table
        .load_preset(UTF8_BORDERS_ONLY)
        .set_content_arrangement(ContentArrangement::Dynamic)
        .set_width(90)
        .set_header(vec![
            Cell::new("Node").add_attribute(Attribute::Bold),
            Cell::new("Type").add_attribute(Attribute::Bold),
            Cell::new("Status").add_attribute(Attribute::Bold),
            Cell::new("Using\ncores").add_attribute(Attribute::Bold),
            Cell::new("Total\ncores").add_attribute(Attribute::Bold),
            Cell::new("User").add_attribute(Attribute::Bold),
            Cell::new("Owner").add_attribute(Attribute::Bold),
        ]);

    let pestat_output = {
        let out = Command::new("sh")
            .args(["-c", "pestat -c | awk '{print $1,$2,$3,$4,$5,$10}'"])
            .output()
            .unwrap()
            .stdout;
        String::from_utf8(out).unwrap()
    };

    for (line, owner) in pestat_output.lines().skip(3).zip(node_owners.iter()) {
        let mut row: Vec<&str> = line.split_whitespace().collect();
        let mut color = Color::Reset;
        if row.len() == 5 {
            row.push("None");
            color = Color::Green;
        }
        if row[5] == std::env::var("USER").unwrap() {
            color = Color::Blue;
        }
        if row[2] == "down*" {
            color = Color::Red;
        }
        row.push(*owner);

        let cells = row.iter().map(|&s| Cell::new(s).fg(color));

        table.add_row(cells);
    }
    println!("{}", table);
}
