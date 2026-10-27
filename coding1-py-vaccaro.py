# %%
import random, time, signal
import matplotlib.pyplot as plt
import threading
import _thread
import numpy as np
import urllib.request
import math
import heapq
from collections import deque, defaultdict, Counter
import csv
import networkx as nx

# ==========================================
# PART 1: PROTEIN SEQUENCES
# ==========================================

# %%
## 1.1. Implement the naive algorithm: for every position of x, compare the window of length k
## against every position of y. State and prove its worst-case bound in n, m and k.
def substring_shared(x, y, n, m, k):
    substring_length = 0
    for i in range(0, n - k + 1):
        for j in range(0, m - k + 1):
            substring_length = 0

            for window in range(k):
                if x[i + window] == y[j + window]:
                    substring_length += 1
                else:
                    break
            if substring_length == k:
                return 1

    return 0


# %%
## 1.2. Implement a hashing-based algorithm: store all k-mers of x in a set, then scan y.
def find_kmer_hashing(x, y, n, m, k):

    if k == 0:
        return "Invalid k value, insert a k of at least two"
    if k > n or k > m:
        return None

    BASE = 131
    MOD = 10**9 + 7

    def val(c):
        return ord(c)

    def first_hash(s, start):
        h = 0
        for t in range(k):
            h = (h * BASE + val(s[start + t])) % MOD
        return h

    power = pow(BASE, k - 1, MOD)

    # Store hashes of all k-mers of x
    x_hashes = {}

    h = first_hash(x, 0)
    x_hashes.setdefault(h, []).append(0)

    for i in range(1, n - k + 1):
        h = (h - val(x[i - 1]) * power) % MOD
        h = (h * BASE + val(x[i + k - 1])) % MOD
        x_hashes.setdefault(h, []).append(i)

    # Scan y
    h = first_hash(y, 0)

    for j in range(m - k + 1):
        if j > 0:
            h = ((h - val(y[j - 1]) * power) * BASE + val(y[j + k - 1])) % MOD
        if h in x_hashes:
            for i in x_hashes[h]:
                if x[i:i + k] == y[j:j + k]:
                    return y[j:j + k]

    return None    


# %%
## 1.3. Benchmark both on random sequences
def time_call(f, args, limit=10):
    """Return elapsed seconds, or None if it exceeds `limit`."""
    def kill():
        _thread.interrupt_main()      # raises KeyboardInterrupt in main thread

    timer = threading.Timer(limit, kill)
    timer.start()
    try:
        start = time.perf_counter()
        f(*args)
        return time.perf_counter() - start
    except KeyboardInterrupt:
        return None
    finally:
        timer.cancel()

def run_benchmark_1_3():
    sizes = [10, 100, 1_000, 10_000, 100_000, 1_000_000]
    ks    = [3, 10, 50]
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    rng = random.Random(0)

    # results[k] = {'naive': [(n, t), ...], 'hash': [(n, t), ...]}
    results = {k: {'naive': [], 'hash': []} for k in ks}

    for n in sizes:
        x = ''.join(rng.choice(alphabet) for _ in range(n))
        y = ''.join(rng.choice(alphabet) for _ in range(n))

        for k in ks:
            if k > n:
                continue

            t_naive = time_call(substring_shared, (x, y, n, n, k))
            t_hash  = time_call(find_kmer_hashing, (x, y, n, n, k))

            results[k]['naive'].append((n, t_naive))
            results[k]['hash'].append((n, t_hash))

            print(f"n={n:>8}  k={k:>3}  "
                  f"naive={t_naive}  hash={t_hash}")

    # ---------- crossover report ----------
    print("\nCrossover points:")
    for k in ks:
        prev_winner = None
        for (n, tn), (_, th) in zip(results[k]['naive'], results[k]['hash']):
            winner = 'hash' if (tn is None or tn > th) else 'naive'
            if prev_winner is not None and winner != prev_winner:
                print(f"  k={k}: crossover near n ≈ {n}  ({prev_winner} -> {winner})")
            prev_winner = winner

    # ---------- plot ----------
    fig, axes = plt.subplots(1, len(ks), figsize=(15, 4.5))

    for ax, k in zip(axes, ks):
        for name, color in [('naive', 'tab:blue'), ('hash', 'tab:orange')]:
            points = [(n, t) for n, t in results[k][name] if t is not None]
            if points:
                ns, ts = zip(*points)
                ax.plot(ns, ts, 'o-', color=color, label=name)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('n = m')
        ax.set_ylabel('time (s)')
        ax.set_title(f'k = {k}')
        ax.legend()
        ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('benchmark.png', dpi=120)
    plt.show()


# ==========================================
# PART 2: DNA SEQUENCES
# ==========================================

# %%
## 2.1. Implement the naive algorithm: for every position and every length l ≥ k, test the substring.
def find_rc_palindromes_naive(t, k):
    n = len(t)
    if k <= 0 or k > n:
        return []
    
    # DNA complement mapping
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    positions = []
    
    # Helper function to check if t[i : i+l] is a reverse-complement palindrome
    def is_rc_pal(i, l):
        # Check characters from outside in
        for j in range(l // 2):
            if comp[t[i + j]] != t[i + l - 1 - j]:
                return False
                
        # For odd lengths, the middle character must be its own complement.
        if l % 2 == 1:
            mid_char = t[i + l // 2]
            if comp[mid_char] != mid_char:
                return False
        return True

    # Iterate over every possible starting position
    for i in range(n - k + 1):
        # Iterate over every possible length >= k
        for l in range(k, n - i + 1):
            if is_rc_pal(i, l):
                positions.append(i)
                break  # We only need to record the position once
                
    return positions


# %%
## 2.2. Implement an algorithm that expands outward from each centre
def find_rc_palindrome_starts(t, k):
    n = len(t)
    if k <= 0 or k > n:
        return []
    
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    starts = set()   # avoid duplicate start positions
    
    # ---------- even-length palindromes ----------
    for c in range(1, n):
        left, right = c - 1, c
        while left >= 0 and right < n:
            if comp[t[left]] != t[right]:
                break
            length = right - left + 1
            if length >= k:
                starts.add(left)
            left -= 1
            right += 1
    
    # ---------- odd-length palindromes ----------
    for c in range(n):
        left, right = c, c
        while left >= 0 and right < n:
            if comp[t[left]] != t[right]:
                break
            length = right - left + 1
            if length >= k:
                starts.add(left)
            left -= 1
            right += 1
    
    return sorted(starts)


# %%
## 2.3. Benchmark both on a real bacterial genome
def run_benchmark_2_3():
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=U00096.3&rettype=fasta&retmode=text"
    print("Downloading E. coli U00096.3 ...")
    fasta = urllib.request.urlopen(url).read().decode()
    genome = "".join(line.strip() for line in fasta.splitlines() if not line.startswith(">"))
    print(f"Genome length: {len(genome):,} bp\n")

    ks = [4, 6, 8, 10]
    results = {}

    for k in ks:
        t0 = time.perf_counter()
        cnt = len(find_rc_palindrome_starts(genome, k))
        elapsed = time.perf_counter() - t0
        results[k] = (cnt, elapsed)
        print(f"k={k:>2}: {cnt:>10,} palindromes  |  {elapsed:6.3f} s")

    print("\n--- Expected vs Observed ---")
    print(f"{'k':>3} {'Observed':>12} {'Expected (random)':>20} {'Ratio (obs/exp)':>18}")
    for k in ks:
        half = math.ceil(k / 2)
        expected = len(genome) * (0.25 ** half)
        obs = results[k][0]
        ratio = obs / expected if expected > 0 else float('inf')
        print(f"{k:>3} {obs:>12,} {expected:>20,.0f} {ratio:>18.2f}")


# ==========================================
# PART 3: GRAPHS
# ==========================================

# %%
## 3.1. Implement BFS and Dijkstra's algorithm
def build_adjlist(G):
    """Return dict: u -> [(v, weight), ...]"""
    adj = {u: [] for u in G.nodes()}
    for u, v, data in G.edges(data=True):
        w = float(data.get("weight", 1.0))
        adj[u].append((v, w))
        if not G.is_directed():
            adj[v].append((u, w))
    return adj

def _reconstruct(parent, v):
    path = []
    while v is not None:
        path.append(v)
        v = parent[v]
    return path[::-1]

def bfs(adj, source):
    """Return (dist, paths): dist[v] = #hops from source; paths[v] = vertex list."""
    dist = {source: 0}
    parent = {source: None}
    q = deque([source])
    while q:
        u = q.popleft()
        for v, _ in adj[u]:
            if v not in dist:
                dist[v] = dist[u] + 1
                parent[v] = u
                q.append(v)
    paths = {v: _reconstruct(parent, v) for v in dist}
    return dist, paths

def dijkstra(adj, source):
    """Return (dist, paths) with minimal total weight from source."""
    dist = {source: 0.0}
    parent = {source: None}
    visited = set()
    pq = [(0.0, source)]
    while pq:
        d, u = heapq.heappop(pq)
        if u in visited:
            continue
        visited.add(u)
        for v, w in adj[u]:
            nd = d + w
            if v not in dist or nd < dist[v]:
                dist[v] = nd
                parent[v] = u
                heapq.heappush(pq, (nd, v))
    paths = {v: _reconstruct(parent, v) for v in dist}
    return dist, paths

def _dists_close(a, b, tol=1e-9):
    if set(a) != set(b):
        return False
    return all(abs(a[k] - b[k]) <= tol for k in a)

def _edge_weight(adj, a, b):
    best = None
    for x, w in adj[a]:
        if x == b and (best is None or w < best):
            best = w
    return best

def _paths_valid(adj, source, dist, paths, weighted, tol=1e-9):
    for v, d in dist.items():
        p = paths[v]
        if not p or p[0] != source or p[-1] != v:
            return False
        total = 0.0
        for a, b in zip(p, p[1:]):
            w = _edge_weight(adj, a, b)
            if w is None:
                return False
            total += w
        if weighted:
            if abs(total - d) > tol:
                return False
        else:
            if len(p) - 1 != d:
                return False
    return True

def run_trial(G, source):
    problems = []
    adj = build_adjlist(G)

    # ---------------- BFS ----------------
    bdist, bpaths = bfs(adj, source)
    nx_bdist = dict(nx.single_source_shortest_path_length(G, source))
    if not _dists_close(bdist, nx_bdist):
        problems.append(("bfs", f"dist mismatch: ours={bdist} nx={nx_bdist}"))
    if not _paths_valid(adj, source, bdist, bpaths, weighted=False):
        problems.append(("bfs", "invalid path"))

    # ---------------- Dijkstra ----------------
    ddist, dpaths = dijkstra(adj, source)
    nx_ddist, _ = nx.single_source_dijkstra(G, source, weight="weight")
    if not _dists_close(ddist, nx_ddist):
        problems.append(("dijkstra", f"dist mismatch: ours={ddist} nx={nx_ddist}"))
    if not _paths_valid(adj, source, ddist, dpaths, weighted=True):
        problems.append(("dijkstra", "invalid path"))

    return problems

# Renamed from 'main' to avoid overwriting the ATP main
def validate_networkx(n_trials=100, seed=12345):
    rng = random.Random(seed)
    total_disagreements = 0
    detail = []

    for t in range(n_trials):
        n = rng.randint(5, 60)
        p = rng.uniform(0.02, 0.85)
        directed = rng.random() < 0.30

        G = nx.gnp_random_graph(
            n, p, seed=rng.randint(0, 2**31 - 1), directed=directed
        )
        for u, v in G.edges():
            G[u][v]["weight"] = rng.uniform(0.1, 10.0)

        source = rng.choice(list(G.nodes()))
        problems = run_trial(G, source)

        if problems:
            total_disagreements += 1
            detail.append((t, n, p, directed, source, problems))

    print(f"Trials: {n_trials}")
    print(f"Graphs with any disagreement: {total_disagreements}")
    if detail:
        for row in detail[:10]:
            print("  ", row)
    else:
        print("No disagreements — BFS and Dijkstra match networkx on all 100 graphs.")


# %%
## 3.2. Benchmarking & Empirical Exponent Estimation
def estimate_exponent(n_vals, times):
    """Fit log(T) = alpha * log(N) + C to estimate exponent alpha."""
    log_n = np.log(n_vals)
    log_t = np.log(times)
    alpha, intercept = np.polyfit(log_n, log_t, 1)
    return alpha

def benchmark_regime(n_values, graph_gen_fn, name):
    bfs_times, dijk_times = [], []

    print(f"\n--- Running Benchmark: {name} ---")
    for n in n_values:
        G = graph_gen_fn(n)
        for u, v in G.edges():
            G[u][v]["weight"] = random.uniform(0.1, 10.0)
        
        adj = build_adjlist(G)
        source = random.choice(list(G.nodes()))

        # Benchmark BFS
        t0 = time.perf_counter()
        _ = bfs(adj, source)
        t_bfs = time.perf_counter() - t0
        bfs_times.append(t_bfs)

        # Benchmark Dijkstra
        t0 = time.perf_counter()
        _ = dijkstra(adj, source)
        t_dijk = time.perf_counter() - t0
        dijk_times.append(t_dijk)

        print(f"n = {n:7d} | BFS: {t_bfs:.5f}s | Dijkstra: {t_dijk:.5f}s")

    alpha_bfs = estimate_exponent(n_values, bfs_times)
    alpha_dijk = estimate_exponent(n_values, dijk_times)

    print(f"\nEmpirical Exponent (alpha in T ~ n^alpha):")
    print(f"  BFS:      alpha = {alpha_bfs:.3f}")
    print(f"  Dijkstra: alpha = {alpha_dijk:.3f}")

    plt.figure(figsize=(7, 5))
    plt.loglog(n_values, bfs_times, 'o-', label=f'BFS (slope ~ {alpha_bfs:.2f})')
    plt.loglog(n_values, dijk_times, 's-', label=f'Dijkstra (slope ~ {alpha_dijk:.2f})')
    plt.xlabel('Number of Vertices (n) [Log Scale]')
    plt.ylabel('Execution Time (seconds) [Log Scale]')
    plt.title(f'Log-Log Benchmark: {name}')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    
def run_benchmarks_3_2():
    n_regime1 = [1000, 5000, 10000, 50000, 100000, 500000, 1000000]
    regime1_gen = lambda n: nx.fast_gnp_random_graph(n, p=4.0/n, seed=42)
    benchmark_regime(n_regime1, regime1_gen, "Regime 1: Mean Degree = 4 (m ~ 2n)")

    n_regime2 = [100, 300, 700, 1500, 3000, 6000, 10000]
    regime2_gen = lambda n: nx.fast_gnp_random_graph(n, p=0.01, seed=42)
    benchmark_regime(n_regime2, regime2_gen, "Regime 2: Fixed Density p = 0.01 (m ~ O(n^2))")

    plt.show()


# %%
## 3.3. ATP-weighted metabolic network
def load_network(path):
    adj = defaultdict(list)
    nodes = set()
    with open(path, newline="", encoding="utf-8-sig") as f:
        rdr = csv.reader(f)
        try:
            next(rdr) # Skip header
        except StopIteration:
            return adj, nodes

        si, ti, wi = 0, 1, 2

        for row in rdr:
            if len(row) <= max(si, ti, wi):
                continue
            u, v = row[si].strip(), row[ti].strip()
            try:
                w = float(row[wi])
            except ValueError:
                continue
            adj[u].append((v, w))
            nodes.update((u, v))
    return adj, nodes

# NOTE: This overwrites the dijkstra from 3.1, as in your original code.
def dijkstra(adj, source):
    dist = {source: 0.0}
    parent = {source: None}
    done = set()
    pq = [(0.0, source)]
    while pq:
        d, u = heapq.heappop(pq)
        if u in done:
            continue
        done.add(u)
        for v, w in adj.get(u, ()):
            nd = d + w
            if v not in dist or nd < dist[v]:
                dist[v] = nd
                parent[v] = u
                heapq.heappush(pq, (nd, v))
    return dist, parent

def bellman_ford(adj, nodes, source):
    INF = float("inf")
    dist = {u: INF for u in nodes}
    parent = {u: None for u in nodes}
    dist[source] = 0.0

    edges = [(u, v, w) for u in nodes for (v, w) in adj.get(u, ())]
    n = len(nodes)

    for _ in range(n - 1):
        changed = False
        for u, v, w in edges:
            if dist[u] != INF and dist[u] + w < dist[v]:
                dist[v] = dist[u] + w
                parent[v] = u
                changed = True
        if not changed:
            break

    neg_cycle = None
    cycle_nodes = set()
    for u, v, w in edges:
        if dist[u] != INF and dist[u] + w < dist[v]:
            x = v
            for _ in range(n):
                x = parent.get(x, x)
            cyc = [x]
            curr = parent.get(x)
            while curr != x and curr is not None:
                cyc.append(curr)
                curr = parent.get(curr)
            if curr == x:
                cyc.append(x)
                cyc.reverse() 
                neg_cycle = cyc
                cycle_nodes.add(v)
                break

    if cycle_nodes:
        for _ in range(n):
            for u, v, w in edges:
                if dist[u] == -INF or dist[u] + w < dist[v]:
                    dist[v] = -INF

    return dist, parent, neg_cycle

def path_to(parent, v):
    p = []
    seen = set()
    while v is not None:
        if v in seen:
            return None  
        seen.add(v)
        p.append(v)
        v = parent.get(v)
    return p[::-1]

def format_path(parent, v):
    p = path_to(parent, v)
    if p is None:
        return "<invalid path: cycle detected>"
    return " -> ".join(p)

def main(path="metabolic_energy_balance_network.csv"):
    t0 = time.perf_counter()
    adj, nodes = load_network(path)
    n_edges = sum(len(v) for v in adj.values())
    print(f"[load] {time.perf_counter()-t0:.4f}s | {len(nodes)} nodes, {n_edges} directed edges")

    n_neg = sum(1 for u in adj for _, w in adj[u] if w < 0)
    n_pos = sum(1 for u in adj for _, w in adj[u] if w > 0)
    n_zero = n_edges - n_neg - n_pos
    print(f"[weights] {n_neg} negative, {n_pos} positive, {n_zero} zero\n")

    src = "Glucose"
    if src not in nodes:
        src = next((u for u in nodes if u.lower() == "glucose"), None)
        if src is None:
            raise SystemExit("Glucose not found in graph")

    targets = ["Pyruvate", "Oxaloacetate", "Serine"]

    print("=== Dijkstra (NOT trustworthy: graph has negative weights) ===")
    dD, pD = dijkstra(adj, src)
    for t in targets:
        if t in dD:
            print(f"  {src} -> {t}: {dD[t]:+g}")
            print(f"    path: {format_path(pD, t)}")
        else:
            print(f"  {src} -> {t}: unreachable")

    print("\n=== Bellman–Ford (correct ground truth) ===")
    t0 = time.perf_counter()
    dB, pB, neg = bellman_ford(adj, nodes, src)
    print(f"  [bellman-ford runtime] {time.perf_counter()-t0:.4f}s")
    for t in targets:
        d_val = dB.get(t, float("inf"))
        if d_val == -float("inf"):
            print(f"  {src} -> {t}: -infinity (reachable from negative cycle)")
        elif d_val < float("inf"):
            print(f"  {src} -> {t}: {d_val:+g}")
            print(f"    path: {format_path(pB, t)}")
        else:
            print(f"  {src} -> {t}: unreachable")

    print("\n=== A concrete place Dijkstra is wrong ===")
    found = False
    for v, d_true in dB.items():
        d_dij = dD.get(v)
        if d_true in (float("inf"), -float("inf")) or d_dij is None:
            continue
        if d_dij - d_true > 1e-9:
            found = True
            print(f"  Vertex: {v}")
            print(f"    Dijkstra says  {d_dij:+g}    path: {format_path(pD, v)}")
            print(f"    True value     {d_true:+g}    path: {format_path(pB, v)}")
            print(
                f"    -> Dijkstra overestimates by {d_dij - d_true:g} because it "
                f"settled {v} early without evaluating subsequent negative edges."
            )
            break
    if not found:
        print("  No discrepancy found for the requested targets/source.")

    print("\n=== Negative cycle reachable from Glucose ===")
    if neg is None:
        print(f"  No negative cycle reachable from {src}.")
    else:
        total = 0.0
        for i in range(len(neg) - 1):
            for uu, w in adj.get(neg[i], ()):
                if uu == neg[i + 1]:
                    total += w
                    break
        print(f"  REACHABLE negative cycle (total weight {total:+g}):")
        print("  " + " -> ".join(neg))


# ==========================================
# MAIN EXECUTION BLOCK (Runs only if executed directly)
# ==========================================
if __name__ == "__main__":
    print("Starting pipeline...")
    
    # 1. Run E. coli Benchmark (Uncomment if you want to run it via terminal)
    # run_benchmark_2_3()
    
    # 2. Run Graph Validation (Uncomment if you want to run it via terminal)
    # validate_networkx()
    
    # 3. Run ATP Network Pipeline
    main()