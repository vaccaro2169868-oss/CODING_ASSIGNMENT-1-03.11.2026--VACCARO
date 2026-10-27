# CODING_ASSIGNMENT-1-03.11.2026--VACCARO

This repository contains the Jupyter Notebook solution for Coding Assignment 1. The assignment is divided into three main algorithmic areas: string matching (protein sequences), DNA palindrome detection, and graph shortest-path algorithms.

## 📋 Overview

The notebook (`CODING_ASSIGNMENT-1(03.11.2026)-VACCARO.ipynb`) contains the following sections:

1. **String Matching (Protein Sequences):** 
   - Compares a naive substring search against a rolling-hash algorithm.
   - Benchmarks both algorithms and identifies the crossover points where hashing becomes faster.
2. **DNA Palindromes:** 
   - Implements naive and center-expansion algorithms to find reverse-complement palindromes (restriction enzyme recognition sites).
   - Benchmarks the algorithms on the real *E. coli* genome (U00096.3) and compares the observed counts against the theoretical random expectation.
3. **Graph Algorithms:** 
   - **3.1 & 3.2:** Implements Breadth-First Search (BFS) and Dijkstra's algorithm from scratch. Validates them against `networkx` on 100 random graphs and benchmarks them across different graph densities.
   - **3.3:** Compares Dijkstra against Bellman-Ford on an ATP-weighted metabolic network. Demonstrates why Dijkstra fails on negative-weight graphs and extracts reachable negative-weight cycles.

## ⚙️ Prerequisites & Installation

To run the code, you will need **Python 3.8+** and Jupyter Notebook (or JupyterLab).

1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/your-repo-name.git
   cd your-repo-name


