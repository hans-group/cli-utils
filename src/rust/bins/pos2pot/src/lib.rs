use std::collections::HashMap;

use inquire::formatter::OptionFormatter;
use inquire::Select;
use serde::{Deserialize, Serialize};

const POTCAR_JSON: &str = include_str!("../data/potcar.json");

#[derive(Serialize, Deserialize)]
pub struct PotcarData {
    pub element: String,
    pub potcar_name: String,
    pub enmax: i32,
    pub recommended: bool,
}

fn get_potcar_list() -> Vec<PotcarData> {
    let potcar_list: Vec<PotcarData> =
        serde_json::from_str(POTCAR_JSON).expect("Error in reading potcar.json");
    potcar_list
}

pub fn get_recommended_potcars(elems: &[String]) -> Vec<String> {
    println!("Getting recommended POTCARs...");
    let potcar_data = get_potcar_list();
    elems
        .iter()
        .map(|e| {
            potcar_data
                .iter()
                .find(|p| &p.element == e && p.recommended)
                .unwrap()
        })
        .map(|p| p.potcar_name.to_string())
        .collect()
}

pub fn get_potcars_from_map(elems: &[String], potcar_map: &HashMap<String, String>) -> Vec<String> {
    println!("Using supplied POTCARs...");
    elems
        .iter()
        .map(|e| potcar_map.get(e).unwrap().to_string())
        .collect()
}

pub fn prompt_potcars(elems: &[String]) -> Vec<String> {
    println!("No POTCARs are entered. Please select POTCARs for each element.");
    let mut selected_potcars: Vec<String> = Vec::new();
    let potcar_data = get_potcar_list();
    let formatter: OptionFormatter<String> = &|i| {
        let words: Vec<String> = i
            .to_string()
            .split_whitespace()
            .map(str::to_string)
            .collect();
        format!("{} selected", words[0])
    };
    for elem in elems {
        let mut potcar_list: Vec<String> = vec![];
        let mut potcar_info_list: Vec<String> = vec![];
        for potcar in &potcar_data {
            if &potcar.element == elem {
                potcar_list.push(potcar.potcar_name.clone());
                let potcar_info: String = format!(
                    "{:<10}{:<15}{:<6}",
                    potcar.potcar_name, potcar.enmax, potcar.recommended
                );
                potcar_info_list.push(potcar_info);
            }
        }
        let potcar_info_table = potcar_info_list
            .iter()
            .zip(potcar_list.iter())
            .collect::<HashMap<_, _>>();

        let potcar_name = Select::new(
            format!(
                "Choose POTCAR for {}\n  {:<10}{:<15}{}",
                elem, "Name", "Cutoff (eV)", "Recommended"
            )
            .as_str(),
            potcar_info_list.clone(),
        )
        .with_formatter(formatter)
        .prompt()
        .unwrap();

        let selected_potcar_name = potcar_info_table.get(&potcar_name).unwrap();
        selected_potcars.push(selected_potcar_name.to_string());
    }
    selected_potcars
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_potcar_recommendation() {
        let elements = vec!["Pt".to_string(), "Cu".to_string()];
        assert_eq!(get_recommended_potcars(&elements), elements);
    }
}
