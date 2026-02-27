/*
This file will contain structures relevant to setting up the proof system.
*/
use ark_ff::Field;
use ark_poly::SparseMultilinearExtension;

pub trait SumCheckProver<F: Field> {
    // Computes the sum so we can have an alleged claim of the functions.
    fn compute_sum(&mut self) -> F;

    // Creates a function that has one variable (meaning it fixes all other variables)
    fn get_verifier_function(&self) -> SparseMultilinearExtension<F>;

    fn fix_variable(&mut self, random_field_element: F);

}

pub trait SumCheckVerifier<F: Field> {
    // Has to check the degree of the function to ensure no one cheats.
    fn verify_degree(&self, fx: &SparseMultilinearExtension<F>) -> bool;

    // Returns a random field element from the verifier
    fn get_random_field_element(&mut self) -> F;

    // Takes as input a multilinear extension and checks that for each field their sum is the claim.
    fn check_claimed_value(&self, fx: &SparseMultilinearExtension<F>) -> bool;

    /// Should ideally take a function by the prover and do all necessary checks
    /// If any fails then it panics, and if everything is good then it returns a random field element.
    fn handle_round(&mut self, fx: &SparseMultilinearExtension<F>) -> F;

    fn set_claim(&mut self, claim: F);

    fn final_check(&self);
}