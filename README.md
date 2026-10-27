# CODING_ASSIGNMENT-1-03.11.2026--VACCARO

Algorithmic Sequence Analysis & Graph Algorithms Benchmark. This repository contains a comprehensive suite of theoretical proofs, algorithmic implementations, and empirical performance benchmarks in Python. The codebase addresses three core computer science and bioinformatics domains: k-mer substring searching, reverse-complement DNA palindrome detection, and shortest-path algorithms on graphs with arbitrary edge weights.

---

## 📋 Table of Contents

- [Overview](#-overview)
- [Module Summary & Algorithms](#-module-summary--algorithms)
  - [1. Substring & k-Mer Search](#1-substring--k-mer-search)
  - [2. Reverse-Complement Palindrome Detection](#2-reverse-complement-palindrome-detection)
  - [3. Graph Shortest Paths & Metabolic Networks](#3-graph-shortest-paths--metabolic-networks)
- [Installation & Requirements](#-installation--requirements)
- [Usage](#-usage)
- [Empirical Benchmarks & Insights](#-empirical-benchmarks--insights)
- [Repository Structure](#-repository-structure)

---

## 🔍 Overview

The notebook provides robust Python implementations alongside mathematical analysis for key algorithmic problems:

* **String Algorithms**: Comparing Naive window search with Polynomial Rolling Hash search across varying sequence lengths ($n$, $m$) and window sizes ($k$).
* **Genomics**: Searching complete bacterial genomes (*E. coli* strain U00096.3) for restriction enzyme recognition sites (reverse-complement palindromes).
* **Graph Theory**: Validating custom BFS and Dijkstra implementations against `networkx` on random graphs ($G(n,p)$) and demonstrating Dijkstra's failure on metabolic networks containing ATP-weighted negative edges and negative-weight cycles.

---

## 🛠️ Module Summary & Algorithms

### 1. Substring & k-Mer Search

* **Problem**: Given two protein sequences $x$ and $y$ over the 20-letter amino acid alphabet of lengths $n$ and $m$, and an integer $k$, return a shared substring of length $k$ if one exists.
* **Implementations**:
  * **Naive Window Comparison** (`substring_shared`):
    * **Worst-Case Time Complexity**: $\mathcal{O}((n - k + 1) \cdot (m - k + 1) \cdot k) = \mathcal{O}(n \cdot m \cdot k)$.
  * **Rolling Hash / Rabin-Karp Search** (`find_kmer_hashing`):
    * **Expected Time Complexity**: $\mathcal{O}(n + m)$ using rolling hash lookup over a hash dictionary.
    * **Worst-Case Time Complexity**: $\mathcal{O}(n \cdot m)$ in the presence of complete hash collisions.

### 2. Reverse-Complement Palindrome Detection

* **Problem**: Locate all starting positions of DNA reverse-complement palindromes of length $l \ge k$ in a sequence $t$ of length $n$.
* **Implementations**:
  * **Naive Substring Testing** (`find_rc_palindromes_naive`): Tests all starting positions $i$ and lengths $l \ge k$. Worst-case complexity: $\mathcal{O}(n^3)$.
  * **Center-Expansion Algorithm** (`find_rc_palindrome_starts`): Expands outward from each potential center position (even and odd lengths) and halts at the first base mismatch.
    * **Time Complexity**: $\mathcal{O}(n + P)$, where $P$ is the total length/number of matched palindromes found.
* **Genomic Validation**: Fetches the *E. coli* U00096.3 genome directly via NCBI Entrez APIs and compares observed palindrome frequencies against random sequence probability expectations ($P_{\text{match}} = 0.25^{\lceil k/2 \rceil}$).

### 3. Graph Shortest Paths & Metabolic Networks

* **Unweighted/Weighted Algorithms**:
  * **Breadth-First Search (BFS)**: Computes single-source unweighted hop distances and shortest paths in $\mathcal{O}(\vert{}V\vert{} + \vert{}E\vert{})$ time.
  * **Dijkstra's Algorithm**: Priority queue-based (`heapq`) implementation for non-negatively weighted graphs running in $\mathcal{O}((\vert{}V\vert{} + \vert{}E\vert{}) \log \vert{}V\vert{})$ time.
  * **Validation**: Verified over 100 random $G(n,p)$ graph trials with zero disagreements compared to NetworkX.
* **Metabolic Network Analysis (`metabolic_energy_balance_network.csv`)**:
  * Evaluates net ATP cost reactions (+1 consumption, -1 production, 0 neutral).
  * **Dijkstra Failure Analysis**: Explains why Dijkstra fails in graphs with negative edge weights due to greedy distance finalization.
  * **Bellman-Ford Algorithm**: Implemented with $\mathcal{O}(\vert{}V\vert{} \cdot \vert{}E\vert{})$ complexity to handle negative edge weights and detect/extract reachable **negative-weight cycles**.

---

## 📦 Installation & Requirements

Ensure you have Python 3.9+ installed along with the following required libraries:

```bash
pip install numpy matplotlib networkx

## 🚀 Usag

Launch Jupyter Notebook or JupyterLab from your terminal:

jupyter notebook "CODING_ASSIGNMENT-1(03.11.2026)-VACCARO.ipynb"


