//! Shared tests for the sum-check provers.
//!
//! Behavioral tests are written once, generic over [`SumCheckProver`], and
//! instantiated for both the naive and the fast prover. Reference sums are
//! computed with plain arkworks operations (adapted from arkworks'
//! `gkr_round_sumcheck` tests).

use ark_bls12_381::Fr;
use ark_ff::{Field, One, Zero};
use ark_poly::univariate::SparsePolynomial;
use ark_poly::MultilinearExtension;
use ark_std::{test_rng, UniformRand};

use crate::sumcheck::{
    FastProver, GKRRound, NaiveProver, StandardVerifier, SumCheckProver, SumCheckVerifier,
};
use crate::util::{index_to_field_element, random_gate};

const DIM: usize = 7;

/// A random sum-check instance together with a random fixed gate label.
/// Deterministic: `test_rng` is seeded, so every call returns the same instance.
fn random_instance() -> (GKRRound<Fr>, Vec<Fr>) {
    let round: GKRRound<Fr> = GKRRound::new_rand(DIM);
    let gate = random_gate(round.gate_labels());
    (round, gate)
}

fn prover<P: SumCheckProver<Fr>>() -> P {
    let (round, gate) = random_instance();
    P::new(round, &gate)
}

/// Both provers, built over the *same* instance and gate.
fn prover_pair() -> (FastProver<Fr>, NaiveProver<Fr>) {
    let (round, gate) = random_instance();
    (FastProver::new(round.clone(), &gate), NaiveProver::new(round, &gate))
}

/// Reference implementation of the initial sum using only arkworks operations.
fn ark_compute_sum<F: Field>(round: &GKRRound<F>, gate: &[F]) -> F {
    let dim = round.vi().num_vars;
    let mult_at_gate = round.mult_predicate().fix_variables(gate);
    let add_at_gate = round.add_predicate().fix_variables(gate);
    let mut sum_xy = F::zero();
    for x in 0..(1 << dim) {
        let f2_x = round.vi()[x];
        let x_index = index_to_field_element(x, dim);
        let mf1_gx = mult_at_gate
            .fix_variables(&x_index)
            .to_dense_multilinear_extension();
        let af1_gx = add_at_gate
            .fix_variables(&x_index)
            .to_dense_multilinear_extension();
        for y in 0..(1 << dim) {
            let add_term = af1_gx[y] * f2_x + af1_gx[y] * round.vj()[y];
            sum_xy += mf1_gx[y] * f2_x * round.vj()[y] + add_term;
        }
    }
    sum_xy
}

/// Reference implementation of the sum after one variable has been fixed to `r`.
fn ark_compute_sum_after_fixing<F: Field>(round: &GKRRound<F>, gate: &[F], r: F) -> F {
    let mult_at_gate = round.mult_predicate().fix_variables(gate);
    let add_at_gate = round.add_predicate().fix_variables(gate);
    let f2 = round.vi();
    let f3 = round.vj();
    let mut sum_xy = F::zero();
    for x in 0..(1 << (f2.num_vars - 1)) {
        let f2_x = f2.fix_variables(&[r])[x];
        let x_index = index_to_field_element(x, f2.num_vars - 1);
        let mf1_gx = mult_at_gate
            .fix_variables(&[r])
            .fix_variables(&x_index)
            .to_dense_multilinear_extension();
        let af1_gx = add_at_gate
            .fix_variables(&[r])
            .fix_variables(&x_index)
            .to_dense_multilinear_extension();
        for y in 0..(1 << f3.num_vars) {
            let add_term = af1_gx[y] * f2_x + af1_gx[y] * f3[y];
            sum_xy += mf1_gx[y] * f2_x * f3[y] + add_term;
        }
    }
    sum_xy
}

fn evaluations_sum(points: &[(usize, Fr)]) -> Fr {
    points.iter().rfold(Fr::zero(), |acc, (_, elem)| *elem + acc)
}

// ---------------------------------------------------------------------------
// Generic behavioral tests, instantiated for both provers.
// ---------------------------------------------------------------------------

fn assert_initial_sum_matches_ark_reference<P: SumCheckProver<Fr>>() {
    let (round, gate) = random_instance();
    let mut prover = P::new(round.clone(), &gate);
    assert_eq!(prover.compute_sum(), ark_compute_sum(&round, &gate));
}

#[test]
fn naive_initial_sum_matches_ark_reference() {
    assert_initial_sum_matches_ark_reference::<NaiveProver<Fr>>();
}

#[test]
fn fast_initial_sum_matches_ark_reference() {
    assert_initial_sum_matches_ark_reference::<FastProver<Fr>>();
}

fn assert_verifier_accepts_initial_round<P: SumCheckProver<Fr>>() {
    let mut prover: P = prover();
    // Degree bound 2 is arbitrary for the simulation; the round polynomials
    // sent here are degree 1.
    let verifier = StandardVerifier::new(2, prover.compute_sum());
    let points = prover.get_verifier_function();
    assert!(verifier.check_claimed_value(&SparsePolynomial::from_coefficients_vec(points)));
}

#[test]
fn naive_verifier_accepts_initial_round() {
    assert_verifier_accepts_initial_round::<NaiveProver<Fr>>();
}

#[test]
fn fast_verifier_accepts_initial_round() {
    assert_verifier_accepts_initial_round::<FastProver<Fr>>();
}

fn assert_verifier_function_sums_to_claim<P: SumCheckProver<Fr>>() {
    let mut prover: P = prover();
    let points = prover.get_verifier_function();
    assert_eq!(prover.compute_sum(), evaluations_sum(&points));
}

#[test]
fn naive_verifier_function_sums_to_claim() {
    assert_verifier_function_sums_to_claim::<NaiveProver<Fr>>();
}

#[test]
fn fast_verifier_function_sums_to_claim() {
    assert_verifier_function_sums_to_claim::<FastProver<Fr>>();
}

/// After fixing the next variable to 0 (resp. 1), the new sum must equal the
/// round polynomial's evaluation g(0) (resp. g(1)).
fn assert_fixing_constant_yields_evaluation<P: SumCheckProver<Fr>>(value: Fr, index: usize) {
    let mut prover: P = prover();
    let verifier_function = prover.get_verifier_function();
    prover.fix_variable(value);
    assert_eq!(prover.compute_sum(), verifier_function[index].1);
}

#[test]
fn naive_fixing_zero_yields_g0() {
    assert_fixing_constant_yields_evaluation::<NaiveProver<Fr>>(Fr::zero(), 0);
}

#[test]
fn fast_fixing_zero_yields_g0() {
    assert_fixing_constant_yields_evaluation::<FastProver<Fr>>(Fr::zero(), 0);
}

#[test]
fn naive_fixing_one_yields_g1() {
    assert_fixing_constant_yields_evaluation::<NaiveProver<Fr>>(Fr::one(), 1);
}

#[test]
fn fast_fixing_one_yields_g1() {
    assert_fixing_constant_yields_evaluation::<FastProver<Fr>>(Fr::one(), 1);
}

fn assert_fixing_matches_ark_reference<P: SumCheckProver<Fr>>() {
    let (round, gate) = random_instance();
    let mut prover = P::new(round.clone(), &gate);

    let r = Fr::rand(&mut test_rng());
    prover.fix_variable(r);

    assert_eq!(prover.compute_sum(), ark_compute_sum_after_fixing(&round, &gate, r));
}

#[test]
fn naive_fixing_matches_ark_reference() {
    assert_fixing_matches_ark_reference::<NaiveProver<Fr>>();
}

#[test]
fn fast_fixing_matches_ark_reference() {
    assert_fixing_matches_ark_reference::<FastProver<Fr>>();
}

// ---------------------------------------------------------------------------
// Agreement between the two provers.
// ---------------------------------------------------------------------------

#[test]
fn provers_agree_on_initial_sum() {
    let (mut fast, mut naive) = prover_pair();
    assert_eq!(naive.compute_sum(), fast.compute_sum());
}

#[test]
fn provers_verifier_functions_agree() {
    let (mut fast, mut naive) = prover_pair();
    let fast_function = fast.get_verifier_function();
    let naive_function = naive.get_verifier_function();

    assert_eq!(evaluations_sum(&naive_function), evaluations_sum(&fast_function));
    assert_eq!(naive_function[0], fast_function[0]);
    assert_eq!(naive_function[1], fast_function[1]);
}

#[test]
fn provers_agree_after_fixing_one() {
    let (mut fast, mut naive) = prover_pair();
    let fast_function = fast.get_verifier_function();
    let naive_function = naive.get_verifier_function();

    fast.fix_variable(Fr::one());
    naive.fix_variable(Fr::one());

    assert_eq!(fast.compute_sum(), fast_function[1].1);
    assert_eq!(naive.compute_sum(), naive_function[1].1);
    assert_eq!(naive_function[1], fast_function[1]);
}

#[test]
fn provers_agree_after_fixing_random() {
    let (mut fast, mut naive) = prover_pair();
    let r = Fr::rand(&mut test_rng());

    fast.fix_variable(r);
    naive.fix_variable(r);

    assert_eq!(fast.compute_sum(), naive.compute_sum());
}

/// Covers both phases: rounds 1..n fix the `x` variables, rounds n+1..2n fix
/// the `y` variables. Regression test for the phase-two initialization, which
/// used to drop the addition term.
#[test]
fn provers_agree_across_all_rounds() {
    let (round, gate) = random_instance();
    let mut fast = FastProver::new(round.clone(), &gate);
    let mut naive = NaiveProver::new(round.clone(), &gate);

    for i in 0..2 * round.gate_labels() {
        let r = Fr::rand(&mut test_rng());
        fast.fix_variable(r);
        naive.fix_variable(r);
        assert_eq!(fast.compute_sum(), naive.compute_sum(), "diverged at round {}", i + 1);
    }
    assert_eq!(fast.fixed_variables().len(), 2 * round.gate_labels());
}

// ---------------------------------------------------------------------------
// Instance construction.
// ---------------------------------------------------------------------------

#[test]
fn random_instance_has_expected_dimensions() {
    let (round, gate) = random_instance();
    assert_eq!(round.mult_predicate().num_vars, 3 * DIM);
    assert_eq!(gate.len(), DIM);
}
