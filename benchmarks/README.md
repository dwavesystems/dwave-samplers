# SA Benchmarks

This folder contains benchmarking scripts for simulated annealing backends.

## Benchmark: runtime vs `num_reads`

Script: `benchmark_sa_reads.py`

Compares only:
- `cpu_sa`
- `fast_cpu_sa`

The script sweeps a list of `num_reads` values, runs both backends on the same
generated Ising problem, and writes:
- CSV summary
- PNG plot

## Usage

From repo root:

```bash
python benchmarks/benchmark_sa_reads.py
```

Custom sweep:

```bash
python benchmarks/benchmark_sa_reads.py \
  --reads "100,1000,5000,10000,20000,50000" \
  --num-variables 100 \
  --edge-probability 0.1 \
  --num-sweeps 1000 \
  --repeats 3
```

## Outputs

Default output paths:
- `benchmarks/results/sa_reads_benchmark.csv`
- `benchmarks/results/sa_reads_benchmark.png`

## Notes

- Requires `matplotlib` to render the plot.
- Uses median runtime in the plot to reduce noise from outliers.
