use regex::Regex;

pub fn read_energies(oszicar: &str) -> Vec<f64> {
    let re_e0 = Regex::new(r"E0=\s+[\d+]?\+?-?.\d+E[\+-]?\d+").unwrap();
    let res = oszicar
        .lines()
        .filter(|&line| line.contains("E0"))
        .map(|line| {
            re_e0
                .captures(line)
                .and_then(|cap| cap.get(0))
                .and_then(|m| m.as_str().split("=").last())
                .and_then(|s| s.trim().parse::<f64>().ok())
        });
    res.map(|x| x.unwrap()).collect()
}
