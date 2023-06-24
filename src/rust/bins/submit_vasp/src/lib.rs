use home::home_dir;
use std::path::Path;
use std::process::Command;

/// Check if a VASP job is properly configured.
/// Checks the following:
/// - POSCAR, INCAR, KPOINTS, POTCAR exists
pub fn check_vasp_missing_files(calcdir: &str) -> Result<(), String> {
    let calcdir = Path::new(calcdir);
    let poscar = calcdir.join("POSCAR");
    let incar = calcdir.join("INCAR");
    let kpoints = calcdir.join("KPOINTS");
    let potcar = calcdir.join("POTCAR");

    let missing_list = vec![poscar, incar, kpoints, potcar]
        .iter()
        .filter(|x| !x.exists())
        .map(|x| x.to_str().unwrap().to_string())
        .collect::<Vec<String>>();

    match missing_list.len() {
        0 => Ok(()),
        _ => Err(format!(
            "The following files are missing: {:?}",
            missing_list
        )),
    }
}

fn get_all_partitions() -> Vec<String> {
    let output = Command::new("sh")
        .arg("-c")
        .arg("sinfo -s | tail -n +2 | awk '{print $1}' | sed 's/*//g'")
        .output()
        .expect("Failed to execute command");
    let output = String::from_utf8(output.stdout).unwrap();
    let partitions = output.trim().split("\n").map(|x| x.to_string()).collect();
    partitions
}

pub fn check_partition(partition: &str) -> Result<String, String> {
    let partitions = get_all_partitions();
    let given_partitions = partition.split(",").collect::<Vec<&str>>();
    match given_partitions
        .iter()
        .all(|x| partitions.contains(&x.to_string()))
    {
        true => Ok(partition.to_string()),
        false => Err(format!(
            "Partition {} does not exist. Available partitions are: {:?}",
            partition, partitions
        )),
    }
}

pub fn get_vasp_bin(vasp_bin_name: &str) -> Result<String, String> {
    let config_dir = home_dir().unwrap().join(".config");
    let config_file = config_dir.join("vasp_bins.json");
    // Check if config file exists
    if !config_file.exists() {
        return Err(
            "Config file '~/.config/vasp_bins.json' does not exist. Please create it.".into(),
        );
    } else {
        let config_file = std::fs::File::open(config_file).map_err(|e| e.to_string())?;
        let config: serde_json::Value =
            serde_json::from_reader(config_file).map_err(|e| e.to_string())?;
        let vasp_bin = config[vasp_bin_name].as_str();
        match vasp_bin {
            Some(vasp_bin) => Ok(vasp_bin.to_string()),
            None => Err(format!(
                "VASP binary '{}' not found in config file.",
                vasp_bin_name
            )),
        }
    }
}
