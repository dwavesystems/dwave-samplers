#!/usr/bin/env python3
"""Benchmark SA runtime vs num_reads for cpu_sa and fast_cpu_sa backends.

Outputs:
  - CSV with raw benchmark summary
  - PNG plot of wall time vs num_reads
"""

from __future__ import annotations

import argparse
import csv
import statistics
import time
from pathlib import Path

import numpy as np
import dimod

from dwave.samplers import SimulatedAnnealingSampler


def build_problem(num_variables: int, edge_probability: float, seed: int) -> dimod.BinaryQuadraticModel:
    rng = np.random.default_rng(seed)

    h = {i: float(rng.uniform(-1.0, 1.0)) for i in range(num_variables)}
    j = {}
    for u in range(num_variables):
        for v in range(u + 1, num_variables):
            if rng.random() < edge_probability:
                j[(u, v)] = float(rng.uniform(-1.0, 1.0))

    return dimod.BinaryQuadraticModel.from_ising(h, j)


def parse_reads(value: str) -> list[int]:
    reads = [int(x.strip()) for x in value.split(",") if x.strip()]
    if not reads:
        raise ValueError("num_reads list cannot be empty")
    if min(reads) <= 0:
        raise ValueError("all num_reads values must be > 0")
    return sorted(reads)


def benchmark(
    sampler: SimulatedAnnealingSampler,
    bqm: dimod.BinaryQuadraticModel,
    backend: str,
    num_reads_values: list[int],
    num_sweeps: int,
    repeats: int,
    warmup_runs: int,
) -> list[dict]:
    # warmup to avoid one-time overhead in measured runs
    for _ in range(warmup_runs):
        sampler.sample(
            bqm,
            num_reads=num_reads_values[0],
            num_sweeps=num_sweeps,
            sa_backend=backend,
            seed=123,
        )

    rows = []
    for num_reads in num_reads_values:
        samples = []
        for i in range(repeats):
            t0 = time.perf_counter()
            sampler.sample(
                bqm,
                num_reads=num_reads,
                num_sweeps=num_sweeps,
                sa_backend=backend,
                seed=1000 + i,
            )
            dt = time.perf_counter() - t0
            samples.append(dt)

        rows.append(
            {
                "backend": backend,
                "num_reads": num_reads,
                "mean_seconds": statistics.mean(samples),
                "median_seconds": statistics.median(samples),
                "min_seconds": min(samples),
                "max_seconds": max(samples),
            }
        )
    return rows


def write_csv(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "backend",
        "num_reads",
        "mean_seconds",
        "median_seconds",
        "min_seconds",
        "max_seconds",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_plot(rows: list[dict], path: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover
        raise RuntimeError(
            "matplotlib is required to generate plots. Install with `pip install matplotlib`."
        ) from exc

    cpu = [r for r in rows if r["backend"] == "cpu_sa"]
    fast = [r for r in rows if r["backend"] == "fast_cpu_sa"]

    cpu = sorted(cpu, key=lambda r: r["num_reads"])
    fast = sorted(fast, key=lambda r: r["num_reads"])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(
        [r["num_reads"] for r in cpu],
        [r["median_seconds"] for r in cpu],
        marker="o",
        label="cpu_sa",
    )
    ax.plot(
        [r["num_reads"] for r in fast],
        [r["median_seconds"] for r in fast],
        marker="o",
        label="fast_cpu_sa",
    )

    ax.set_title("Simulated Annealing Runtime vs num_reads")
    ax.set_xlabel("num_reads")
    ax.set_ylabel("median runtime (seconds)")
    ax.grid(True, alpha=0.3)
    ax.legend()

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reads",
        default="100,500,1000,5000,10000,20000,50000",
        help="Comma-separated num_reads values",
    )
    parser.add_argument("--num-variables", type=int, default=100)
    parser.add_argument("--edge-probability", type=float, default=0.1)
    parser.add_argument("--num-sweeps", type=int, default=1000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmup-runs", type=int, default=1)
    parser.add_argument("--problem-seed", type=int, default=7)
    parser.add_argument(
        "--csv-out",
        default="benchmarks/results/sa_reads_benchmark.csv",
        help="Output CSV path",
    )
    parser.add_argument(
        "--plot-out",
        default="benchmarks/results/sa_reads_benchmark.png",
        help="Output plot path",
    )
    args = parser.parse_args()

    num_reads_values = parse_reads(args.reads)
    sampler = SimulatedAnnealingSampler()
    bqm = build_problem(args.num_variables, args.edge_probability, args.problem_seed)

    rows = []
    for backend in ("cpu_sa", "fast_cpu_sa"):
        rows.extend(
            benchmark(
                sampler=sampler,
                bqm=bqm,
                backend=backend,
                num_reads_values=num_reads_values,
                num_sweeps=args.num_sweeps,
                repeats=args.repeats,
                warmup_runs=args.warmup_runs,
            )
        )

    csv_path = Path(args.csv_out)
    plot_path = Path(args.plot_out)
    write_csv(rows, csv_path)
    write_plot(rows, plot_path)

    print(f"Wrote CSV:  {csv_path}")
    print(f"Wrote plot: {plot_path}")


if __name__ == "__main__":
    main()
