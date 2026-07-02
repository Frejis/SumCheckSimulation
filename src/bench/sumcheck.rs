//! Standalone sum-check benchmark: naive vs. fast prover on random instances
//! of increasing dimension. Results are written to `sum-check.xlsx`.

use ark_bls12_381::Fr;
use ark_ff::Field;
use ark_poly::univariate::SparsePolynomial;
use rust_xlsxwriter::{ColNum, RowNum, Workbook};

use crate::sumcheck::{
    FastProver, GKRRound, NaiveProver, StandardVerifier, SumCheckProver, SumCheckVerifier,
};
use crate::timing::{timed, Track};
use crate::util::random_gate;

const OUTPUT_FILE: &str = "sum-check.xlsx";
const MAX_VARIABLES: usize = 10;

/// Trials per dimension as used for the report benchmarks: dimensions 0-5 are
/// skipped, 6-10 run 300 trials and everything above runs 30.
fn trials_per_dimension() -> Vec<usize> {
    [vec![0; 6], vec![300; 5], vec![30; 10]].concat()
}

/// Runs the sum-check benchmark configuration used in the report and saves
/// the results to `sum-check.xlsx`. Note: takes roughly 5 hours in total.
pub fn run_sum_check_benchmark() {
    simulate_and_save(MAX_VARIABLES, &trials_per_dimension());
}

fn simulate_and_save(max_variables: usize, trials: &[usize]) {
    let mut naive_results: Vec<Vec<Track>> = Vec::new();
    let mut fast_results: Vec<Vec<Track>> = Vec::new();

    for (dimension, &trial_count) in trials.iter().enumerate().take(max_variables) {
        let mut naive_trials = Vec::new();
        let mut fast_trials = Vec::new();
        for _ in 0..trial_count {
            let naive_time = simulate_instance::<NaiveProver<Fr>, Fr>(dimension);
            println!("Finished running naive for dimension {dimension}");
            println!("Time taken for naive prover: {:?}", naive_time.prover());
            naive_trials.push(naive_time);

            let fast_time = simulate_instance::<FastProver<Fr>, Fr>(dimension);
            println!("Finished running fast for dimension {dimension}");
            println!("Time taken for fast prover: {:?}", fast_time.prover());
            fast_trials.push(fast_time);
        }
        naive_results.push(naive_trials);
        fast_results.push(fast_trials);
    }

    save_results(&fast_results, &naive_results, max_variables);
}

/// Simulates one full sum-check execution over a random instance of the given
/// dimension and returns the accumulated prover/verifier times.
fn simulate_instance<P: SumCheckProver<F>, F: Field>(dimension: usize) -> Track {
    let round: GKRRound<F> = GKRRound::new_rand(dimension);
    let gate = random_gate(round.gate_labels());
    let mut track = Track::new();

    let (mut prover, elapsed) = timed(|| P::new(round.clone(), &gate));
    track.add_prover_time(elapsed);
    let (claim, elapsed) = timed(|| prover.compute_sum());
    track.add_prover_time(elapsed);

    // Degree bound 2 is arbitrary for the simulation; the round polynomials
    // sent here are degree 1.
    let (mut verifier, elapsed) = timed(|| StandardVerifier::new(2, claim));
    track.add_verifier_time(elapsed);

    let mut verifier_function = Vec::new();
    for _ in 0..2 * round.vj.num_vars {
        let (points, elapsed) = timed(|| prover.get_verifier_function());
        track.add_prover_time(elapsed);
        verifier_function = points.clone();

        let (r_i, elapsed) =
            timed(|| verifier.handle_round(&SparsePolynomial::from_coefficients_vec(points)));
        track.add_verifier_time(elapsed);

        let (_, elapsed) = timed(|| prover.fix_variable(r_i));
        track.add_prover_time(elapsed);

        let (claim, elapsed) = timed(|| prover.compute_sum());
        track.add_prover_time(elapsed);

        let (_, elapsed) = timed(|| verifier.set_claim(claim));
        track.add_verifier_time(elapsed);
    }

    let (_, elapsed) = timed(|| {
        verifier.final_check(
            &gate,
            round.add_predicate(),
            round.mult_predicate(),
            round.vi.clone(),
            verifier_function,
        )
    });
    track.add_verifier_time(elapsed);

    track
}

fn save_results(fast: &[Vec<Track>], naive: &[Vec<Track>], max_variables: usize) {
    let mut workbook = Workbook::new();
    let worksheet = workbook.add_worksheet();

    for i in 0..max_variables {
        let col = (i * 4) as ColNum;
        worksheet.write(0, col, format!("Fast Prover, {i} variables")).unwrap();
        worksheet.write(0, col + 1, format!("Fast Verifier, {i} variables")).unwrap();
        worksheet.write(0, col + 2, format!("Naive Prover, {i} variables")).unwrap();
        worksheet.write(0, col + 3, format!("Naive Verifier, {i} variables")).unwrap();
    }

    for (i, (fast_trials, naive_trials)) in fast.iter().zip(naive).enumerate() {
        for (j, (fast_trial, naive_trial)) in fast_trials.iter().zip(naive_trials).enumerate() {
            let row = (1 + j) as RowNum;
            let col = (i * 4) as ColNum;
            worksheet.write(row, col, fast_trial.prover().as_secs_f64()).unwrap();
            worksheet.write(row, col + 1, fast_trial.verifier().as_secs_f64()).unwrap();
            worksheet.write(row, col + 2, naive_trial.prover().as_secs_f64()).unwrap();
            worksheet.write(row, col + 3, naive_trial.verifier().as_secs_f64()).unwrap();
        }
    }

    workbook.save(OUTPUT_FILE).unwrap();
    println!("Saved data to {OUTPUT_FILE}");
}
