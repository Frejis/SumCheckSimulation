# SumCheckSimulation
This is a simulation of the sum-check protocol. It will naturally follow my own implementation of the naive implementation compared with the prefix suffix sum-check protocol, both of which are compared to the sum-check protocol designed in arkworks in the Rust ecosystem.

# Plan
1. Implement the naive sumcheck protocol for a GKR round. Look at arkworks for inspiration on parameter values and setup.
    1. The idea is to compare the two, so my interface should follow theirs to allow for easy comparisons.
1. Implement the prefix-sufix sumcheck protocol for a GKR round. Again look at arkworks to ensure a matching interface such that benchmarking is smooth.
1. Automatic figures generated that illustrates the benchmark so it is easy to compare.
1. Create a makefile such that you can run make benchmark which will generate the benchmarks given some prerequisites (like having rust installed etc.).
