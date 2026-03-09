use std::time::Duration;
use ark_bls12_381::Fr;
use ark_ff::Field;
use ark_std::test_rng;
use structures::data_structures::SumCheckProver;
use crate::gkr::gkr_protocol::simulate_gkr_circuit;
use crate::gkr::gkr_round::GKRRound;
use crate::provers::{fast, naive};
use crate::provers::fast::FastProver;
use crate::provers::libra::Libra;
use crate::provers::naive::NaiveProver;
use crate::structures::circuit_structures::GKRCircuit;

mod util;
pub mod provers;
pub mod structures;
pub mod verifiers;
pub mod gkr;
pub mod benchmarking;

fn main() {
    simulate_circuit_and_print_results();
}

fn simulate_circuit_and_print_results() {
    let config = BenchmarkConfig::default();
    let data_path = benchmarking::circuit_data::default_benchmark_data_path(&config.layer_sizes, config.trials);

    // Try to load pre-generated circuits, or generate them if they don't exist
    let circuit_set = load_or_generate_circuits(&config, &data_path);

    // Benchmark each prover on the same circuit set
    let fast_result = benchmark_prover_on_circuit_set("Fast", &config, &circuit_set, fast_prover_ctor);
    let naive_result = benchmark_prover_on_circuit_set("Naive", &config, &circuit_set, naive_prover_ctor);
    let libra_result = benchmark_prover_on_circuit_set("Libra", &config, &circuit_set, libra_prover_ctor);

    print_benchmark_report(&config, &[fast_result, naive_result, libra_result]);
}

#[derive(Clone)]
struct BenchmarkConfig {
    layer_sizes: Vec<usize>,
    trials: usize,
}

impl Default for BenchmarkConfig {
    fn default() -> Self {
        Self {
            layer_sizes: vec![1, 2, 4, 8, 16, 32, 64, 128],
            trials: 1,
        }
    }
}

struct BenchmarkResult {
    backend: &'static str,
    total_prover_time: Duration,
    total_verifier_time: Duration,
    trials: usize,
}

impl BenchmarkResult {
    fn average_prover_time(&self) -> Duration {
        average_duration(self.total_prover_time, self.trials)
    }

    fn average_verifier_time(&self) -> Duration {
        average_duration(self.total_verifier_time, self.trials)
    }

    fn total_protocol_time(&self) -> Duration {
        self.total_prover_time + self.total_verifier_time
    }

    fn average_protocol_time(&self) -> Duration {
        average_duration(self.total_protocol_time(), self.trials)
    }
}

fn average_duration(total: Duration, samples: usize) -> Duration {
    if samples == 0 {
        Duration::ZERO
    } else {
        total / samples as u32
    }
}

fn fast_prover_ctor(round: GKRRound<Fr>, random_gate: Vec<Fr>) -> FastProver<Fr> {
    FastProver::new(round, &random_gate)
}

fn naive_prover_ctor(round: GKRRound<Fr>, random_gate: Vec<Fr>) -> NaiveProver<Fr> {
    NaiveProver::new(round, &random_gate)
}

fn libra_prover_ctor(round: GKRRound<Fr>, random_gate: Vec<Fr>) -> Libra<Fr> {
    Libra::new(&round, random_gate)
}

/// Load pre-generated circuits from disk, or generate and save them if they don't exist
fn load_or_generate_circuits(
    config: &BenchmarkConfig,
    data_path: &std::path::Path,
) -> benchmarking::circuit_data::BenchmarkCircuitSet {
    use benchmarking::circuit_data::BenchmarkCircuitSet;

    if data_path.exists() {
        println!("Loading pre-generated circuits from {:?}...", data_path);
        match BenchmarkCircuitSet::load_from_file(data_path) {
            Ok(set) => {
                println!("Successfully loaded {} circuits", set.len());
                return set;
            }
            Err(e) => {
                println!("Warning: Failed to load circuits: {}", e);
                println!("Regenerating circuits...");
            }
        }
    } else {
        println!("No pre-generated circuits found at {:?}", data_path);
        println!("Generating {} new circuits...", config.trials);
    }

    generate_and_save_circuits(config, data_path)
}

/// Generate new random circuits and save them to disk
fn generate_and_save_circuits(
    config: &BenchmarkConfig,
    data_path: &std::path::Path,
) -> benchmarking::circuit_data::BenchmarkCircuitSet {
    use benchmarking::circuit_data::BenchmarkCircuitSet;

    let mut rng = test_rng();
    let mut circuit_set = BenchmarkCircuitSet::new();

    println!("Generating circuits with layer sizes: {:?}", config.layer_sizes);
    for i in 0..config.trials {
        let circuit: GKRCircuit<Fr> = GKRCircuit::random(config.layer_sizes.as_slice(), &mut rng);
        circuit_set.add_circuit(&circuit);
        
        if (i + 1) % 10 == 0 || i + 1 == config.trials {
            println!("  Generated {}/{} circuits", i + 1, config.trials);
        }
    }

    println!("Saving circuits to {:?}...", data_path);
    match circuit_set.save_to_file(data_path) {
        Ok(_) => println!("Successfully saved circuits"),
        Err(e) => println!("Warning: Failed to save circuits: {}", e),
    }

    circuit_set
}

/// Benchmark a specific prover on a pre-generated set of circuits
fn benchmark_prover_on_circuit_set<P>(
    backend_name: &'static str,
    config: &BenchmarkConfig,
    circuit_set: &benchmarking::circuit_data::BenchmarkCircuitSet,
    prover_ctor: fn(GKRRound<Fr>, Vec<Fr>) -> P,
) -> BenchmarkResult
where
    P: SumCheckProver<Fr>,
{
    use std::time::Instant;

    println!("\nBenchmarking {} prover...", backend_name);
    
    let mut total_prover_time = Duration::ZERO;
    let mut total_verifier_time = Duration::ZERO;

    for i in 0..config.trials {
        let circuit: GKRCircuit<Fr> = circuit_set
            .get_circuit(i)
            .expect("Circuit index out of bounds");
        
        let start = Instant::now();
        let (prover_time, verifier_time) = simulate_gkr_circuit(circuit, prover_ctor);
        let elapsed = start.elapsed();
        
        total_prover_time += prover_time;
        total_verifier_time += verifier_time;
        
        if (i + 1) % 10 == 0 || i + 1 == config.trials {
            println!("  Completed {}/{} trials (last trial: {:?})", i + 1, config.trials, elapsed);
        }
    }

    BenchmarkResult {
        backend: backend_name,
        total_prover_time,
        total_verifier_time,
        trials: config.trials,
    }
}

fn benchmark_backend<P>(
    backend: &'static str,
    config: &BenchmarkConfig,
    prover_ctor: fn(GKRRound<Fr>, Vec<Fr>) -> P,
) -> BenchmarkResult
where
    P: SumCheckProver<Fr>,
{
    let mut rng = test_rng();
    let mut total_prover_time = Duration::ZERO;
    let mut total_verifier_time = Duration::ZERO;

    for _ in 0..config.trials {
        let circuit = GKRCircuit::random(config.layer_sizes.as_slice(), &mut rng);
        let (prover_time, verifier_time) = simulate_gkr_circuit(circuit, prover_ctor);
        total_prover_time += prover_time;
        total_verifier_time += verifier_time;
    }

    BenchmarkResult {
        backend,
        total_prover_time,
        total_verifier_time,
        trials: config.trials,
    }
}

fn print_benchmark_report(config: &BenchmarkConfig, results: &[BenchmarkResult]) {
    println!("================ GKR Circuit Benchmark ================");
    println!("Layers: {:?}", config.layer_sizes);
    println!("Trials per backend: {}", config.trials);

    for result in results {
        println!("-------------------------------------------------------");
        println!("Backend: {}", result.backend);
        println!("Total prover time: {:?}", result.total_prover_time);
        println!("Total verifier time: {:?}", result.total_verifier_time);
        println!("Total protocol time: {:?}", result.total_protocol_time());
        println!("Avg prover time/trial: {:?}", result.average_prover_time());
        println!("Avg verifier time/trial: {:?}", result.average_verifier_time());
        println!("Avg protocol time/trial: {:?}", result.average_protocol_time());
    }

    println!("=======================================================");
}

#[cfg(test)]
mod generic_tests {
    use std::time::{Duration, Instant};
    use ark_bls12_381::Fr;
    use ark_ff::Field;
    use ark_poly::MultilinearExtension;
    use crate::structures::data_structures::{SumCheckProver, SumCheckVerifier};
    use crate::fast::FastProver;
    use crate::gkr::gkr_round::GKRRound;
    use crate::naive::NaiveProver;
    use crate::util::random_gate;
    use crate::verifiers::standard_verifier::StandardVerifier;

    #[test]
    fn test_verifier_first_round() {
        // 7 variables 3 seconds.... 8 variables 23!!! seconds!??! :OOO whaaa

        let gkr_round: GKRRound<Fr> = GKRRound::new_rand();
        let random_gate = random_gate(gkr_round.gate_labes());
        let mut prover = NaiveProver::new(gkr_round.clone(), &random_gate);
        // Now we test the g_func gives what we expect
        let verifier = StandardVerifier::new(3, prover.compute_sum(), gkr_round);
        let verifier_func = prover.get_verifier_function();

        assert!(verifier.check_claimed_value(&verifier_func));
    }

    #[test]
    fn simulate_two_rounds_naive() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand_var_size(5);
        let random_gate = random_gate(gkr_round.gate_labes());
        let prover = NaiveProver::new(gkr_round.clone(), &random_gate);
        _compare_verifier_sum(gkr_round, prover);
    }

    #[test]
    fn simulate_two_rounds_fast() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand_var_size(5);
        let random_gate = random_gate(gkr_round.gate_labes());
        let prover = FastProver::new(gkr_round.clone(), &random_gate);
        _compare_verifier_sum(gkr_round, prover);
    }

    fn _compare_verifier_sum<F: Field, P: SumCheckProver<F>>(gkr_round: GKRRound<F>, mut prover: P) {
        let mut verifier = StandardVerifier::new(3, prover.compute_sum(), gkr_round.clone());
        let mut prover_time_spent = Duration::ZERO;
        let mut verifier_time_spent = Duration::ZERO;
        for _ in 0..gkr_round.gate_labes() {

            let verifier_func = prover.get_verifier_function();
            let time = Instant::now();
            assert!(verifier.check_claimed_value(&verifier_func));
            let time_diff = time.elapsed();
            verifier_time_spent += time_diff;

            let rand_var = verifier.get_random_field_element();
            prover.fix_variable(rand_var);

            let time = Instant::now();
            let new_claim = prover.compute_sum();
            let time_diff = time.elapsed();
            prover_time_spent += time_diff;
            verifier.set_claim(new_claim);
        }
        println!("Prover time: {:?}", prover_time_spent);
        println!("Verifier time: {:?}", verifier_time_spent);
    }

    #[test]
    fn test_generated_round() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand_var_size(8);
        assert_eq!(gkr_round.vi.num_vars, 8);
        assert_eq!(gkr_round.vj.num_vars, 8);
        assert_eq!(gkr_round.mult().num_vars(), 24);
    }
}