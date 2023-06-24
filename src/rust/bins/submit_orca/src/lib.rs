use std::process::Command;

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
