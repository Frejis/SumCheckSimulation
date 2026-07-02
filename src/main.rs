use std::env;
use std::process::ExitCode;

use sum_check_simulation::bench::{run_gkr_benchmark, run_sum_check_benchmark};

const USAGE: &str = "usage: sum-check-simulation [sumcheck | gkr [--small]]

  sumcheck       run the sum-check benchmark (writes sum-check.xlsx, ~5h)
  gkr            run the GKR benchmark on the large random circuit
                 (writes gkr_circuit.xlsx)
  gkr --small    run the GKR benchmark on the small figure circuit";

fn main() -> ExitCode {
    let args: Vec<String> = env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        None | Some("sumcheck") => run_sum_check_benchmark(),
        Some("gkr") => run_gkr_benchmark(args.iter().any(|arg| arg == "--small")),
        _ => {
            eprintln!("{USAGE}");
            return ExitCode::FAILURE;
        }
    }
    ExitCode::SUCCESS
}
