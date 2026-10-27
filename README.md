# CODING ASSIGNMENT-1 (03.11.2026)-VACCARO

This repository contains the Python Code solution for Coding Assignment 1. Also, it contains a Jupyter Notebook that allows output and plot visualization. The assignment is divided into three main algorithmic areas: string matching (protein sequences), DNA palindrome detection, and graph shortest-path algorithms.

## 📋 Overview

The Python File with the algorithms(`coding1-py-vaccaro`) contains the following sections:

1. **String Matching (Protein Sequences):** 
   - Compares a naive substring search against a rolling-hash algorithm.
   - Benchmarks both algorithms and identifies the crossover points where hashing becomes faster.
2. **DNA Palindromes:** 
   - Implements naive and center-expansion algorithms to find reverse-complement palindromes (restriction enzyme recognition sites).
   - Benchmarks the algorithms on the real *E. coli* genome (U00096.3) and compares the observed counts against the theoretical random expectation.
3. **Graph Algorithms:** 
   - **3.1 & 3.2:** Implements Breadth-First Search (BFS) and Dijkstra's algorithm from scratch. Validates them against `networkx` on 100 random graphs and benchmarks them across different graph densities.
   - **3.3:** Compares Dijkstra against Bellman-Ford on an ATP-weighted metabolic network (`metabolic_energy_balance_network.csv`). Demonstrates why Dijkstra fails on negative-weight graphs and extracts reachable negative-weight cycles.
  
The Jupyter Notebook File (`coding1-ipynb-vaccaro`) contains the same code in the Python one, but allows the user to visualize the plots and the results of the analysis.
NOTE: the Python file doesn't print the output because the tests run is easier and the algortihms do not have to be runned again.

The Python File to run tests (`tests1-ipynb-vaccaro`) imports the functions from the other py file and performs the tests provided.

## 📂 Repository Structure

```python
├── coding1-py-vaccaro.py              # Main algorithms (Parts 1, 2, 3)
├── coding1-ipynb-vaccaro.ipynb        # Same code, with plot/output visualization
├── tests1-ipynb-vaccaro.py            # Unit tests importing the main module
├── metabolic_energy_balance_network.csv  # ATP-weighted metabolic graph (Part 3.3)
└── README.md
```

## 🛠️ Prerequisites

- Python ≥ 3.8

- pip (for installing third-party packages)

- Internet connection (only required to download the E. coli U00096.3 genome during run_benchmark_2_3())

**Standard library used (no installation required)**
random, time, signal, threading, _thread, math, heapq, urllib.request, csv, importlib, unittest — all bundled with Python.

**Third-party libraries**
The main module requires the following external libraries:

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
The test file additionally imports:

import importlib
import unittest
import networkx as nx

## ⚙️ Installation

Clone the repository:

bash
git clone https://github.com/<your-username>/<your-repo>.git
cd <your-repo>
Create a virtual environment (recommended):

bash
python -m venv .venv

# Linux / macOS
source .venv/bin/activate

# Windows (PowerShell)
.venv\Scripts\Activate.ps1
Install the dependencies:

bash
pip install --upgrade pip
pip install matplotlib numpy networkx
Or, if you prefer a one-liner:

bash
pip install matplotlib numpy networkx jupyter
jupyter is only needed if you plan to run the notebook version.

Verify the installation:

python
python -c "import matplotlib, numpy, networkx; print('All good!')"

## ▶️ Usage

Run the main pipeline (Part 3.3 by default)
bash
python coding1-py-vaccaro.py
This loads metabolic_energy_balance_network.csv, computes the shortest paths from Glucose to Pyruvate / Oxaloacetate / Serine using both Dijkstra and Bellman-Ford, prints where Dijkstra fails, and extracts a reachable negative-weight cycle.

Enable additional sections
Uncomment the relevant lines inside the if __name__ == "__main__": block of coding1-py-vaccaro.py:

python
if __name__ == "__main__":
    # run_benchmark_2_3()   # E. coli U00096.3 palindrome benchmark
    # validate_networkx()   # 100 random-graph validation vs networkx
    main()                  # ATP-weighted metabolic network
Run the tests
bash
python tests1-ipynb-vaccaro.py
or, if you prefer the unittest discovery style:

bash
python -m unittest tests1-ipynb-vaccaro
Explore the results visually
bash
jupyter notebook coding1-ipynb-vaccaro.ipynb
All plots (log-log benchmarks, crossover points, exponent fits) are generated inside the notebook.


