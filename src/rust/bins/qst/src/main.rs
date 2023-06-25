use qst::{Job, State};

use clap::Parser;
use comfy_table::modifiers::UTF8_SOLID_INNER_BORDERS;
use comfy_table::presets::UTF8_FULL;
use comfy_table::{Attribute, Cell, ContentArrangement, Table};
use std::cmp::Ordering;
use std::str::FromStr;

fn get_jobs(user: &str) -> Vec<Job> {
    let jobinfo: String = {
        let out = std::process::Command::new("sh")
            .args(["-c", "scontrol show job"])
            .output()
            .unwrap()
            .stdout;
        String::from_utf8(out).unwrap()
    };

    jobinfo
        .split("\n\n")
        .filter(|&x| !x.is_empty())
        .map(|x| Job::from_str(x).unwrap())
        .filter(|x| x.username == user)
        .collect()
}

fn print_table(user: Option<&str>, num_print: Option<i64>, status: &str) {
    let user_from_env = std::env::var("USER").expect("No USER found.");
    let user = match user {
        Some(user) => user,
        None => user_from_env.as_str(),
    };
    // Filter jobs
    let (n_jobs, jobs): (usize, Vec<Job>) = {
        let all_jobs = get_jobs(user);
        let jobs = match State::from_str(status) {
            Ok(state) => all_jobs
                .into_iter()
                .filter(|job| job.state == state)
                .collect(),
            Err(_) => all_jobs,
        };

        let n_jobs = jobs.len();
        let num_print = num_print.unwrap_or(n_jobs as i64);
        let jobs = match num_print.cmp(&0) {
            Ordering::Equal => jobs.into_iter().collect(),
            Ordering::Less => jobs
                .into_iter()
                .skip((n_jobs as i64 + num_print) as usize)
                .collect(),
            Ordering::Greater => jobs.into_iter().take(num_print as usize).collect(),
        };
        (n_jobs, jobs)
    };

    let mut table = Table::new();
    table
        .load_preset(UTF8_FULL)
        .apply_modifier(UTF8_SOLID_INNER_BORDERS)
        .set_content_arrangement(ContentArrangement::Dynamic)
        .set_width(110)
        .set_header(vec![
            Cell::new("Job ID").add_attribute(Attribute::Bold),
            Cell::new("Job name").add_attribute(Attribute::Bold),
            Cell::new("State").add_attribute(Attribute::Bold),
            Cell::new("Partition").add_attribute(Attribute::Bold),
            Cell::new("Num.\nnodes").add_attribute(Attribute::Bold),
            Cell::new("Num.\ntasks").add_attribute(Attribute::Bold),
            Cell::new("Elapsed\ntime").add_attribute(Attribute::Bold),
        ]);

    for job in jobs {
        if job.username == user {
            let color = job.state.color();
            let row = vec![
                Cell::new(job.id).fg(color),
                Cell::new(job.jobname).fg(color),
                Cell::new(job.state.to_str()).fg(color),
                Cell::new(job.partition).fg(color),
                Cell::new(job.numnodes).fg(color),
                Cell::new(job.numtasks).fg(color),
                Cell::new(job.runtime).fg(color),
            ];

            table.add_row(row);
        }
    }

    println!("{}", table);
    println!("Total {} jobs", n_jobs);
}

/// Print status of jobs on the HPC cluster.
#[derive(Parser, Debug)]
#[command(author,version,about,long_about = None)]
struct Args {
    /// The user whose jobs to be printed. Defaults to the current user.
    #[arg(short, long, default_value = None)]
    user: Option<String>,
    /// Number of jobs to be printed. Defaults to all jobs.
    #[arg(short, long, default_value = None, allow_hyphen_values = true)]
    num_print: Option<i64>,
    /// Select status
    #[arg(short, long, default_value = "all")]
    status: String,
}

fn main() {
    let args = Args::parse();
    print_table(args.user.as_deref(), args.num_print, &args.status);
}
