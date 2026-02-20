/*
This file will contain structures relevant to setting up the proof system.
*/
use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, SparseMultilinearExtension};

pub trait Prover<F: Field> {
    // Computes the sum so we can have an alleged claim of the functions.
    fn compute_sum(&self) -> F;

    // Creates a function that has one variable (meaning it fixes all other variables)
    fn get_verifier_function(&self) -> DenseMultilinearExtension<F>;

    fn fix_variable(&mut self, random_field_element: F);

}

pub trait Verifier<F: Field> {
    // Has to check the degree of the function to ensure no one cheats.
    fn verify_degree(fx: DenseMultilinearExtension<F>);

    // Returns a random field element from the verifier
    fn get_random_field_element() -> F;

    // Takes as input a multilinear extension and checks that for each field their sum is the claim.
    fn check_claimed_value(gx: DenseMultilinearExtension<F>);
}

pub trait GKRRoundProver<F: Field> {
    /*
    mult: SparseMultilinearExtension<F>, // mle of mult for round k with (r, i , j)
    vi: DenseMultilinearExtension<F>, // mle of v_(k-1)(i)
    vj: DenseMultilinearExtension<F>, // mle of v_(k-1)(i)
    r: [F], // The gate "r" that is fixed.
    */
    fn set_mult(func: SparseMultilinearExtension<F>);
    fn get_mult() -> SparseMultilinearExtension<F>;
    fn set_vi(func: DenseMultilinearExtension<F>);
    fn get_vi() -> DenseMultilinearExtension<F>;
    fn set_vj(func: DenseMultilinearExtension<F>);
    fn get_vj() -> DenseMultilinearExtension<F>;

    fn set_gate(gate: [F]);
    fn get_gate() -> [F];

}