use clap::Parser;
use serde::Serialize;
use std::{
    error::Error,
    io::{self, Write},
    process::exit,
};
use submit_orca::check_partition;
use tinytemplate::TinyTemplate;

const SLURM_SCRIPT_TEMPLATE_HEAD: &str = include_str!("../data/submit_orca.sh.head.template");
const SLURM_SCRIPT_BODY: &str = include_str!("../data/submit_orca.sh.body");

#[derive(Debug, Serialize)]
struct Context {
    job_name: String,
    partition: String,
    num_tasks_per_node: usize,
    node: String,
    use_nodelist: bool,
    remove_tdir: bool,
}

#[derive(Parser)]
#[command(name = "submit_orca", author = "Minjoon Hong")]
struct Args {
    // Calculation directory.
    #[arg(short = 'd', long, default_value = ".")]
    calc_dir: String,
    /// Job name.
    #[arg(short = 'J', long)]
    job_name: Option<String>,
    /// Partition name.
    #[arg(short, long, default_value = "g4")]
    partition: String,
    /// Number of tasks per node.
    #[arg(short, long, default_value = "32")]
    num_tasks_per_node: usize,
    /// Name of node. Optional.
    #[arg(long)]
    node: Option<String>,
    /// Remove temporary directory.
    #[arg(short, long)]
    remove_tdir: bool,
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
    template.add_template("slurm_script", SLURM_SCRIPT_TEMPLATE_HEAD)?;

    let args = Args::parse();
    if std::path::Path::new(&args.calc_dir).exists() {
        std::env::set_current_dir(&args.calc_dir)?;
    } else {
        eprintln!("Directory {} does not exist.", args.calc_dir);
        exit(1);
    }
    let partition = check_partition(&args.partition)?;
    let (use_nodelist, node) = match args.node {
        Some(node) => (true, node),
        None => (false, String::from("")),
    };
    let job_name = match args.job_name {
        Some(job_name) => job_name,
        None => get_curr_dirname(),
    };
    let context = Context {
        job_name,
        partition,
        num_tasks_per_node: args.num_tasks_per_node,
        node,
        use_nodelist,
        remove_tdir: args.remove_tdir,
    };
    let rendered = template.render("slurm_script", &context)?;
    let mut slurm_script = std::fs::File::create("submit_orca.sh")?;
    slurm_script.write_all(rendered.as_bytes())?;
    slurm_script.write_all(SLURM_SCRIPT_BODY.as_bytes())?;

    if !args.dry {
        let output = std::process::Command::new("sbatch")
            .arg("submit_orca.sh")
            .output()?;
        io::stdout().write_all(&output.stdout).unwrap();
        io::stderr().write_all(&output.stderr).unwrap();
        exit(output.status.code().unwrap_or(0));
    }
    Ok(())
}
