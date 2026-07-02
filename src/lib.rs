//! Simulation and benchmarking of the sum-check protocol and the GKR protocol.
//!
//! The crate is split into four layers:
//! - [`sumcheck`]: the sum-check protocol itself — the prover/verifier traits,
//!   a naive prover, a linear-time prover and the standard verifier.
//! - [`gkr`]: layered arithmetic circuits and the GKR protocol (prover,
//!   verifier and a driver that runs one sum-check per circuit layer).
//! - [`timing`]: small instrumentation types used to time prover and verifier
//!   work separately.
//! - [`bench`]: the benchmark harnesses that produce the `.xlsx` result files.

pub mod bench;
pub mod gkr;
pub mod sumcheck;
pub mod timing;
pub mod util;
