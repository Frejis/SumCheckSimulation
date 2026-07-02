//! Instrumentation used by the simulations to time prover and verifier work
//! separately.

use std::time::{Duration, Instant};

/// Runs `f` and returns its result together with the elapsed wall-clock time.
pub fn timed<T>(f: impl FnOnce() -> T) -> (T, Duration) {
    let start = Instant::now();
    let result = f();
    (result, start.elapsed())
}

/// Accumulated prover/verifier wall-clock time for one protocol run (or one
/// layer of a GKR run).
#[derive(Clone, Default)]
pub struct Track {
    prover: Duration,
    verifier: Duration,
}

impl Track {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn new_times(prover: Duration, verifier: Duration) -> Self {
        Self { prover, verifier }
    }

    pub fn add_prover_time(&mut self, time: Duration) {
        self.prover += time;
    }

    pub fn add_verifier_time(&mut self, time: Duration) {
        self.verifier += time;
    }

    pub fn prover(&self) -> Duration {
        self.prover
    }

    pub fn verifier(&self) -> Duration {
        self.verifier
    }
}

/// Prover/verifier times per GKR layer for one full protocol run.
/// The final input-layer check is appended as an extra entry.
pub struct AnalysisResult {
    time_per_layer: Vec<Track>,
}

impl AnalysisResult {
    pub fn new() -> Self {
        Self { time_per_layer: Vec::new() }
    }

    pub fn add_time_per_layer(&mut self, track: Track) {
        self.time_per_layer.push(track);
    }

    pub fn get_time_for_layer(&self, layer: usize) -> &Track {
        &self.time_per_layer[layer]
    }
}

impl Default for AnalysisResult {
    fn default() -> Self {
        Self::new()
    }
}
