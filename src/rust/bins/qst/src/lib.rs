use comfy_table::Color;
use lazy_static::lazy_static;
use regex::Regex;
use std::error::Error;
use std::str::FromStr;

lazy_static! {
    static ref ID_REGEX: Regex = Regex::new(r"JobId=\d+").unwrap();
    static ref JOBNAME_REGEX: Regex = Regex::new(r"JobName=[/\w\-_\d.]+").unwrap();
    static ref USERNAME_REGEX: Regex = Regex::new(r"UserId=\w+").unwrap();
    static ref STATE_REGEX: Regex = Regex::new(r"JobState=\w+").unwrap();
    static ref PARTITION_REGEX: Regex = Regex::new(r"Partition=\w+").unwrap();
    static ref NUMNODES_REGEX: Regex = Regex::new(r"NumNodes=\d+").unwrap();
    static ref NUMTASKS_REGEX: Regex = Regex::new(r"NumCPUs=\d+").unwrap();
    static ref RUNTIME_REGEX: Regex = Regex::new(r"RunTime=[\-\d]+:\d+:\d+").unwrap();
    static ref WORKDIR_REGEX: Regex = Regex::new(r"WorkDir=[/\w\-_\d.]+").unwrap();
}

#[derive(Debug)]
pub struct Job {
    pub id: String,
    pub jobname: String,
    pub username: String,
    pub state: State,
    pub partition: String,
    pub numnodes: String,
    pub numtasks: String,
    pub runtime: String,
    pub workdir: String,
}

impl FromStr for Job {
    type Err = Box<dyn Error>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let id = Job::extract_value(s, &ID_REGEX);
        let jobname = Job::extract_value(s, &JOBNAME_REGEX);
        let username = Job::extract_value(s, &USERNAME_REGEX);
        let state = Job::extract_value(s, &STATE_REGEX);
        let partition = Job::extract_value(s, &PARTITION_REGEX);
        let numnodes = Job::extract_value(s, &NUMNODES_REGEX);
        let numtasks = Job::extract_value(s, &NUMTASKS_REGEX);
        let runtime = Job::extract_value(s, &RUNTIME_REGEX);
        let workdir = Job::extract_value(s, &WORKDIR_REGEX);

        Ok(Job {
            id,
            jobname,
            username,
            state: State::from_str(&state)?,
            partition,
            numnodes,
            numtasks,
            runtime,
            workdir,
        })
    }
}

impl Job {
    fn extract_value(s: &str, re: &Regex) -> String {
        let caps = re.captures(s).unwrap();
        let val = &caps[0].split('=').last();
        match val {
            Some(s) => s.to_string(),
            None => String::from("None"),
        }
    }
}

#[derive(Debug, PartialEq)]
pub enum State {
    Running,
    Completed,
    Pending,
    Failed,
    Cancelled,
}

impl FromStr for State {
    type Err = Box<dyn Error>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "RUNNING" | "R" => Ok(State::Running),
            "COMPLETED" | "C" => Ok(State::Completed),
            "PENDING" | "P" | "Q" => Ok(State::Pending),
            "FAILED" | "F" => Ok(State::Failed),
            "CANCELLED" | "c" => Ok(State::Cancelled),
            // Return error if state is not valid
            _ => Err("Invalid state".into()),
        }
    }
}

impl State {
    pub fn to_str(&self) -> &str {
        match self {
            State::Running => "RUNNING",
            State::Completed => "COMPLETED",
            State::Pending => "PENDING",
            State::Failed => "FAILED",
            State::Cancelled => "CANCELLED",
        }
    }
    pub fn color(&self) -> Color {
        match self {
            State::Running => Color::Grey,
            State::Completed => Color::Green,
            State::Pending => Color::Reset,
            State::Failed => Color::Red,
            State::Cancelled => Color::Yellow,
        }
    }
}
