//! GKR benchmark: runs the full protocol with the naive and the fast
//! sum-check prover over the same circuit and writes the per-layer times to
//! `gkr_circuit.xlsx`.

use ark_bls12_381::Fr;
use ark_std::test_rng;
use rust_xlsxwriter::{ColNum, RowNum, Workbook};

use crate::gkr::{
    compute_predicates, GKRCircuit, GKRDriver, GKRProver, GKRVerifier, InputLayer,
};
use crate::sumcheck::{FastProver, NaiveProver, SumCheckProver};
use crate::timing::AnalysisResult;

const OUTPUT_FILE: &str = "gkr_circuit.xlsx";
const TRIALS: usize = 30;

/// Benchmarks the GKR protocol. `small` selects the small hand-built figure
/// circuit instead of the large random one used in the report.
pub fn run_gkr_benchmark(small: bool) {
    let (layer_sizes, circuit, input) = if small { figure_circuit() } else { random_circuit() };

    let mut naive_results: Vec<AnalysisResult> = Vec::new();
    let mut fast_results: Vec<AnalysisResult> = Vec::new();
    for trial in 0..TRIALS {
        println!("Running trial {}/{TRIALS}", trial + 1);
        naive_results.push(simulate::<NaiveProver<Fr>>(&circuit, &input));
        fast_results.push(simulate::<FastProver<Fr>>(&circuit, &input));
    }

    save_results(&fast_results, &naive_results, &layer_sizes);
}

/// Runs one full GKR protocol execution with the given sum-check prover.
fn simulate<P: SumCheckProver<Fr>>(
    circuit: &GKRCircuit<Fr>,
    input: &InputLayer<Fr>,
) -> AnalysisResult {
    let predicates = compute_predicates(circuit, input);
    let prover = GKRProver::new(circuit, input, predicates.clone());
    let verifier = GKRVerifier::new(input.clone(), predicates);
    let mut driver = GKRDriver::new(prover, verifier, circuit.clone(), input.clone());
    driver.run_circuit::<P>()
}

/// The large random circuit benchmarked in the report.
fn random_circuit() -> (Vec<usize>, GKRCircuit<Fr>, InputLayer<Fr>) {
    let layer_sizes = vec![1, 1024, 1024, 1024, 1024];
    let circuit = GKRCircuit::random(&layer_sizes, &mut test_rng());
    let input = InputLayer::random(*layer_sizes.last().unwrap());
    (layer_sizes, circuit, input)
}

/// The small hand-built circuit from the report's figure.
fn figure_circuit() -> (Vec<usize>, GKRCircuit<Fr>, InputLayer<Fr>) {
    let layer_sizes = vec![2, 4, 4];
    let circuit = GKRCircuit::figure_circuit();
    let input = InputLayer::random(*layer_sizes.last().unwrap());
    (layer_sizes, circuit, input)
}

fn save_results(fast: &[AnalysisResult], naive: &[AnalysisResult], layer_sizes: &[usize]) {
    let mut workbook = Workbook::new();
    let worksheet = workbook.add_worksheet();

    for i in 0..layer_sizes.len() {
        let col = (i * 5) as ColNum;
        worksheet.write(0, col, "Variables/Layer").unwrap();
        worksheet.write(0, col + 1, format!("Fast Prover layer {i}")).unwrap();
        worksheet.write(0, col + 2, format!("Fast Verifier layer {i}")).unwrap();
        worksheet.write(0, col + 3, format!("Naive Prover layer {i}")).unwrap();
        worksheet.write(0, col + 4, format!("Naive Verifier layer {i}")).unwrap();
    }

    for (i, layer_size) in layer_sizes.iter().enumerate() {
        worksheet.write((i + 1) as RowNum, 0, layer_size.ilog2()).unwrap();
        for (j, (fast_run, naive_run)) in fast.iter().zip(naive).enumerate() {
            let row = (j + 1) as RowNum;
            let col = (i * 5) as ColNum;
            let fast_layer = fast_run.get_time_for_layer(i);
            let naive_layer = naive_run.get_time_for_layer(i);
            worksheet.write(row, col + 1, fast_layer.prover().as_secs_f64()).unwrap();
            worksheet.write(row, col + 2, fast_layer.verifier().as_secs_f64()).unwrap();
            worksheet.write(row, col + 3, naive_layer.prover().as_secs_f64()).unwrap();
            worksheet.write(row, col + 4, naive_layer.verifier().as_secs_f64()).unwrap();
        }
    }

    workbook.save(OUTPUT_FILE).unwrap();
    println!("Saved data to {OUTPUT_FILE}");
}
