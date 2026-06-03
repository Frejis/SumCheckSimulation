# SumCheckSimulation
This is a simulation of the sum-check protocol. 
It will naturally follow my own implementation of the naive implementation compared with the prefix suffix sum-check protocol.

The project also contains a simulation of the GKR protocol.
The project is still in the development phase and as such running it can be quite quirky.

If you wish to run a Sum-check simualtion go into the main.rs file in the src folder. and overwrite the main function to:
```rust 
fn main() {
    run_sum_check_config();
}
```
This will run the Sum-check simulation for the config used when benchmarking in the report. Note on my machine it took
roughly 5 hours in total.

To run the GKR Simulation overwrite the main function to
```rust
fn main() {
    benchmark_gkr();
}
```
Running this will benchmark the large circuit in the program. If you wish to simulate the small circuit go to the function
`benchmark_gkr` and overwrite this:
```rust
fn benchmark_gkr() {
    let (layers, random_circuit, input_layer) = random_circuit();
    ...
}
```
with:
```rust
fn benchmark_gkr() {
    let (layers, random_circuit, input_layer) = figure_circuit();
    ...
}
```

# Details on files created.
Running the Sum-check will create a `sum-check.xlsx` file that contains the results.
Running the GKR simulation will create a file named `gkr_circuit.xlsx` containing the results.