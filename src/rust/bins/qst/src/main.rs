use qst::Job;

use clap::Parser;
use comfy_table::modifiers::UTF8_SOLID_INNER_BORDERS;
use comfy_table::presets::UTF8_FULL;
use comfy_table::*;
use std::str::FromStr;

fn print_table(user: Option<&str>, num_print: Option<usize>) {
    let user_from_env = std::env::var("USER").expect("No USER found.");
    let user = match user {
        Some(user) => user,
        None => user_from_env.as_str(),
    };
    let jobinfo = {
        let out = std::process::Command::new("sh")
            .args(["-c", "scontrol show job"])
            .output()
            .unwrap()
            .stdout;
        String::from_utf8(out).unwrap()
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

    let mut rows = vec![];
    for l in jobinfo.split("\n\n").filter(|&x| !x.is_empty()) {
        let job: Job = Job::from_str(l).unwrap();
        if job.username == user {
            let color = match job.state.as_str() {
                "RUNNING" => Color::Grey,
                "COMPLETED" => Color::Green,
                _ => Color::Reset,
            };

            let row = vec![
                Cell::new(job.id).fg(color),
                Cell::new(job.jobname).fg(color),
                Cell::new(job.state).fg(color),
                Cell::new(job.partition).fg(color),
                Cell::new(job.numnodes).fg(color),
                Cell::new(job.numtasks).fg(color),
                Cell::new(job.runtime).fg(color),
            ];
            rows.push(row);
        }
    }
    // Print only latest `num_print` jobs.
    let n_row = rows.len();
    for row in rows.into_iter().rev().take(num_print.unwrap_or(n_row)) {
        table.add_row(row);
    }
    println!("{}", table);
}

/// Print status of jobs on the HPC cluster.
#[derive(Parser, Debug)]
#[command(author,version,about,long_about = None)]
struct Args {
    /// The user whose jobs to be printed. Defaults to the current user.
    #[arg(short, long, default_value = None)]
    user: Option<String>,
    /// Number of jobs to be printed. Defaults to all jobs.
    #[arg(short, long, default_value = None)]
    num_print: Option<usize>,
}

fn main() {
    let args = Args::parse();
    print_table(args.user.as_deref(), args.num_print);
}
