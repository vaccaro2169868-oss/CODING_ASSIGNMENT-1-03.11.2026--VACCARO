# CODING ASSIGNMENT-1 (03.11.2026) — VACCARO

This repository contains the Python code solution for **Coding Assignment 1**. It also contains a Jupyter Notebook that lets you visualise the plots and outputs. The assignment is split into three algorithmic areas: **string matching (protein sequences)**, **DNA palindrome detection**, and **graph shortest-path algorithms**.

## 📋 Overview

The main Python file (`coding1-py-vaccaro.py`) is organised as follows:

1. **String Matching (Protein Sequences)**
   - Naive substring search vs. a rolling-hash algorithm (`substring_shared`, `find_kmer_hashing`).
   - Benchmarks both and reports the crossover point at which hashing becomes faster (`run_benchmark_1_3`).
2. **DNA Palindromes**
   - Naive and centre-expansion algorithms for reverse-complement palindromes / restriction-enzyme recognition sites (`find_rc_palindromes_naive`, `find_rc_palindrome_starts`).
   - Benchmarks both on the real *E. coli* genome (U00096.3) and compares observed counts to the theoretical random expectation (`run_benchmark_2_3`).
3. **Graph Algorithms**
   - **3.1 & 3.2** — BFS and Dijkstra implemented from scratch, validated against `networkx` on 100 random graphs (`validate_networkx`) and benchmarked across two density regimes (`run_benchmarks_3_2`).
   - **3.3** — Dijkstra vs. Bellman-Ford on the ATP-weighted metabolic network (`metabolic_energy_balance_network.csv`); demonstrates Dijkstra's failure on negative weights and extracts reachable negative-weight cycles (`main`).

The notebook (`coding1-ipynb-vaccaro.ipynb`) contains the same code but runs every section so you can see plots and printed results interactively.

> **Note on execution** — Every section is implemented as a standalone function. Only `main()` (Part 3.3) is invoked from the `__main__` block, so running the script directly executes the ATP-network pipeline only. The other benchmarks must be invoked explicitly — either by calling them from `main()`, or by importing them in a small driver script or notebook cell.

The test file (`tests1-ipynb-vaccaro.py`) imports the algorithms from the main module and runs the provided unit tests.

## 📂 Repository structure

```
.
├── coding1-py-vaccaro.py                 # Main algorithms (Parts 1, 2, 3)
├── coding1-ipynb-vaccaro.ipynb           # Same code with plot/output visualisation
├── tests1-ipynb-vaccaro.py               # Unit tests importing the main module
├── metabolic_energy_balance_network.csv  # ATP-weighted metabolic graph (Part 3.3)
└── README.md
```

## 🛠️ Prerequisites

- Python ≥ 3.8
- `pip`
- Internet connection — **only** required by `run_benchmark_2_3()` to download the *E. coli* U00096.3 genome from NCBI.

### Standard library

Imports in `coding1-py-vaccaro.py`:

```
random, time, threading, _thread,
math, heapq, urllib.request, csv,
collections (deque, defaultdict, Counter)
```

Additional imports in `tests1-ipynb-vaccaro.py`:

```
importlib, unittest
```

### Third-party libraries

`coding1-py-vaccaro.py`:

```python
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
```

`tests1-ipynb-vaccaro.py`:

```python
import networkx as nx
```

## ⚙️ Installation

1. **Clone the repository**

   ```bash
   git clone https://github.com/vaccaro2169868-oss/CODING_ASSIGNMENT-1-03.11.2026--VACCARO.git
   cd CODING_ASSIGNMENT-1-03.11.2026--VACCARO
   ```

2. **Create and activate a virtual environment** (recommended)

   ```bash
   python -m venv .venv
   # Linux / macOS
   source .venv/bin/activate
   # Windows (PowerShell)
   .venv\Scripts\Activate.ps1
   ```

3. **Install dependencies**

   ```bash
   pip install --upgrade pip
   pip install matplotlib numpy networkx
   ```

   Add `jupyter` if you want to open the notebook:

   ```bash
   pip install matplotlib numpy networkx jupyter
   ```

4. **Verify**

   ```bash
   python -c "import matplotlib, numpy, networkx; print('All good!')"
   ```

## ▶️ Usage

### Run the default pipeline 

```bash
python coding1-py-vaccaro.py
```

### Run the benchmarks explicitly & explore interactively

Every other section is implemented as a function but not invoked from `__main__`. To run them, create a small driver script (or use a notebook cell) that imports and calls them. To make it easier you can safely run the Notebook to explicit results and plots

```bash
jupyter coding1-ipynb-vaccaro.ipynb
```

### Run the tests

```bash
python tests1-ipynb-vaccaro.py
```

## 📊 Outputs produced

| Artefact | Produced by |
|----------|-------------|
| `benchmark.png` | `run_benchmark_1_3()` (Part 1.3) |
| Console table of observed vs. expected palindromes | `run_benchmark_2_3()` (Part 2.3) |
| Validation report over 100 random graphs | `validate_networkx()` (Part 3.1) |
| Log-log plots with fitted exponent α | `run_benchmarks_3_2()` (Part 3.2) |
| Dijkstra-vs-Bellman-Ford report + negative cycle | `main()` (Part 3.3) |



- **`dijkstra` is redefined in Part 3.3.** This is intentional: the Part-3.3 version works with the `defaultdict(list)` adjacency produced by `load_network`, while the Part-3.1 version expects the `build_adjlist` dict form. The Part-3.1 version is still used by `run_trial` / `validate_networkx`, since those call `build_adjlist` first.

- **`time_call` timeout.** Uses `threading.Timer` + `_thread.interrupt_main`, which is CPython-specific and not sub-second accurate. A call exceeding the 10 s limit returns `None` and is filtered out of the plots.

- **Reproducibility.** All random seeds are fixed — `random.Random(0)` (Part 1.3), `seed=42` (Part 3.2), `seed=12345` (Part 3.1) — so benchmarks and validation reports reproduce across machines.

- **Offline use.** `run_benchmark_2_3()` requires internet access to fetch U00096.3. If you're offline, skip it or cache the FASTA locally and adapt the loader.
