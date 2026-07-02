//! Benchmark harnesses. Each writes its results to an `.xlsx` file in the
//! working directory.

mod gkr;
mod sumcheck;

pub use gkr::run_gkr_benchmark;
pub use sumcheck::run_sum_check_benchmark;
