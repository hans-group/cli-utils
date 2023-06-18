use clap::Parser;
use serde::Serialize;
use std::{
    error::Error,
    io::{self, Write},
    process::exit,
};
use submit_vasp::{check_partition, check_vasp_missing_files, get_vasp_bin};
use tinytemplate::TinyTemplate;

const SLURM_SCRIPT_TEMPLATE: &str = include_str!("../data/submit_vasp.sh.template");

#[derive(Debug, Serialize)]
struct Context {
    job_name: String,
    partition: String,
    num_nodes: usize,
    num_tasks_per_node: usize,
    nodelist: String,
    use_nodelist: bool,
    exclude: String,
    use_exclude: bool,
    vasp_bin: String,
}

#[derive(Parser)]
#[command(name = "submit_orca", author = "Minjoon Hong")]
struct Args {
    /// Calculation directory.
    #[arg(short = 'd', long, default_value = ".")]
    calc_dir: String,
    /// Job name.
    #[arg(short = 'J', long)]
    job_name: Option<String>,
    /// Partition name.
    #[arg(short, long, default_value = "g1")]
    partition: String,
    /// Number of nodes.
    #[arg(short = 'N', long, default_value = "1")]
    num_nodes: usize,
    /// Number of tasks per node.
    #[arg(short, long, default_value = "16")]
    num_tasks_per_node: usize,
    /// Name of node. Optional.
    #[arg(short = 'w', long)]
    nodelist: Option<String>,
    /// List of nodes to exclude
    #[arg(long)]
    exclude: Option<String>,
    /// VASP binary name. Default is 6.3.2/std
    #[arg(long, default_value = "6.3.2/std")]
    vasp_bin: String,
    /// Do not submit job
    #[arg(long)]
    dry: bool,
}

fn get_curr_dirname() -> String {
    let curr_dir = std::env::current_dir().unwrap();
    let curr_dir = curr_dir.to_str().unwrap();
    let curr_dir = curr_dir.split('/').last().unwrap();
    String::from(curr_dir)
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut template = TinyTemplate::new();
    template.add_template("slurm_script", SLURM_SCRIPT_TEMPLATE)?;

    let args = Args::parse();
    if std::path::Path::new(&args.calc_dir).exists() {
        std::env::set_current_dir(&args.calc_dir)?;
    } else {
        eprintln!("Directory {} does not exist.", args.calc_dir);
        exit(1);
    }
    check_vasp_missing_files(".")?;

    let partition = check_partition(&args.partition)?;
    let (use_nodelist, nodelist) = match args.nodelist {
        Some(node) => (true, node),
        None => (false, String::from("")),
    };
    let (use_exclude, exclude) = match args.exclude {
        Some(exclude) => (true, exclude),
        None => (false, String::from("")),
    };

    match (use_nodelist, use_exclude) {
        (true, true) => {
            eprintln!("Cannot use both nodelist and exclude.");
            exit(1);
        }
        _ => (),
    }

    let job_name = match args.job_name {
        Some(job_name) => job_name,
        None => get_curr_dirname(),
    };
    let vasp_bin = get_vasp_bin(&args.vasp_bin)?;
    let context = Context {
        job_name,
        partition,
        num_nodes: args.num_nodes,
        num_tasks_per_node: args.num_tasks_per_node,
        nodelist,
        use_nodelist,
        exclude,
        use_exclude,
        vasp_bin,
    };
    let rendered = template.render("slurm_script", &context)?;
    let mut slurm_script = std::fs::File::create("submit_vasp.sh")?;
    slurm_script.write_all(rendered.as_bytes())?;

    if !args.dry {
        let output = std::process::Command::new("sbatch")
            .arg("submit_vasp.sh")
            .output()?;
        io::stdout().write_all(&output.stdout).unwrap();
        io::stderr().write_all(&output.stderr).unwrap();
        exit(output.status.code().unwrap_or(0));
    }
    Ok(())
}
