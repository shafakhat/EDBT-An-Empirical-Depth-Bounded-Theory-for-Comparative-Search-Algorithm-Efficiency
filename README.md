# EDBT:UnifiedTreeExperiment

A Java implementation of the **EDBT (Extended Dynamic Behavioural Tree)** framework — an accuracy-aware extension for benchmarking and comparing tree and graph search algorithm efficiency. The framework empirically validates theoretical complexity bounds across 10 classical search algorithms using randomised grid environments.

> **Author:** Shafakhatullah Khan Mohammed  
> **Version:** 1.0 | **Since:** 2025-10-18

---

## Overview

`UnifiedTreeExperimentL` runs a comprehensive experimental suite that:

1. Executes **10 search algorithms** across multiple tree depths and random seeds
2. Computes **EDBT metrics** — complexity profile bounds, pruning accuracy, convergence operators, and efficiency indices
3. Validates theoretical bounds (ACP, PCCT) empirically across all algorithm–depth combinations
4. Renders an **interactive Swing GUI** showing node expansions (log scale) and ABEI bar charts

---

## Algorithms Benchmarked

| Algorithm | Type | Notes |
|---|---|---|
| **BFS** | Uninformed | Breadth-first search |
| **Dijkstra** | Uninformed | Optimal single-source shortest path |
| **A\*** | Informed | Weight = 1.0 (optimal) |
| **WA\*(1.5)** | Informed | Weighted A*, ε-suboptimal |
| **WA\*(2.0)** | Informed | Weighted A*, higher suboptimality |
| **GBFS** | Informed | Greedy best-first search |
| **BiDijkstra** | Bidirectional | Bidirectional Dijkstra |
| **RBFS** | Memory-bounded | Recursive best-first search |
| **IDA\*** | Iterative deepening | Iterative deepening A* |
| **AlphaBeta** | Game tree | Alpha-Beta pruning on synthetic tree |

---

## EDBT Framework Components

### 1. Complexity Profile
Estimates the empirical branching exponent **β** via log-linear regression and derives three complexity bounds:

| Symbol | Meaning |
|---|---|
| γ⁺ | Upper bound on node expansions |
| γ⁻ | Lower bound on node expansions |
| γᵘ | Obstacle-aware complexity bound |

### 2. Accuracy Function (Ψ-based)
Computes two pruning accuracy measures at each expansion step using Ψ-values (f-scores) of all successors:

| Metric | Description |
|---|---|
| **Δˢ** (deltaS) | Solution retention ratio — how well the best successor is preserved |
| **Δᶜ** (deltaC) | Pruning ratio — fraction of successors discarded |

### 3. Convergence Operator
Given β and Δᵈ, derives:
- **d\*** — minimum depth for ε-convergence
- **q\*** — quality of solution at convergence

### 4. ACP Bound (Accuracy-Complexity Product)
Formalises the tradeoff between node expansions and solution accuracy:
> **ACP = N̄ · Δᵈ ≥ γ⁻(β̂, d)**  

Verified at every depth for every algorithm.

### 5. PCCT Certificate (Pruning-Certified Convergence Threshold)
Derives a theoretical lower bound on solution quality for each algorithm class (exact, WA\*, GBFS), verified against empirical Δˢ values.

### 6. Efficiency Indices

| Index | Formula | Description |
|---|---|---|
| **BEI** | (U · q) / (N · log_b N) | Branching Efficiency Index |
| **ABEI** | (U · q · Δᵈ) / (N · log_b N) | Accuracy-weighted BEI |

---

## Configuration

| Parameter | Default | Description |
|---|---|---|
| `NUM_TRIALS` | 30 | Independent trials per depth |
| `DEPTHS` | 4, 6, 8, 10, 12 | Search depths swept |
| `BRANCHING` | 4 | Branching factor |
| `OBSTACLE_RATE` | 0.15 | Fraction of grid cells blocked |
| `AB_MAX_DEPTH` | 9 | Max depth for Alpha-Beta tree |
| `EPSILON` | 0.05 | Convergence tolerance ε |

---

## Requirements

- Java 11 or later
- No external dependencies — uses Java standard library and Swing (included in JDK)

---

## Build & Run

```bash
# Compile
javac UnifiedTreeExperimentL.java

# Run
java UnifiedTreeExperimentL
```

On launch the program:
1. Prints full experiment tables to the console
2. Opens a **1400×700 Swing window** with two panels:
   - Left: Node expansions vs depth (log scale) for all 10 algorithms
   - Right: ABEI bar chart per algorithm and depth

---

## Console Output Tables

| Table | Contents |
|---|---|
| **Main sweep** | Per algorithm per depth: N̄ ± std, Δᵈ, δs, BEI, ABEI, ACP vs γ⁻, q\* |
| **β Parameters** | Estimated β̂, α̂, R² for each algorithm |
| **Complexity Profile** | γ⁻, γᵘ, N̄, γ⁺ per depth, in-bounds check |
| **ACP Bound** | N̄ · Δᵈ ≥ γ⁻ verification (✓/✗) per algorithm × depth |
| **PCCT Table** | Δˢ vs certified lower bound per algorithm × depth |
| **ABEI vs BEI** | Side-by-side comparison of both efficiency indices |
| **Convergence** | d\*, q\* for each algorithm |
| **BEI Table** | BEI summary across all depths |
| **Held-out MAPE** | Mean Absolute Percentage Error of β-based predictions |
| **Ranking (d=12)** | Algorithms ranked by ABEI at maximum depth |

---

## Project Structure

```
UnifiedTreeExperimentL.java
│
├── AlgResult               # Per-algorithm result record (all EDBT metrics)
├── EDBTMetrics             # Per-run metrics container
├── Grid                    # Randomised weighted grid environment
│
├── complexityProfile()     # Computes γ⁻, γᵘ, γ⁺ from β and depth
├── accuracyMeasures()      # Computes Δˢ and Δᶜ from Ψ-value sets
├── convergenceOperator()   # Derives d* and q* from β and Δᵈ
├── computePathDeltaD()     # Path-level Δᵈ product along solution path
├── pcctCertificate()       # PCCT lower bound per algorithm class
├── computeABEI()           # Accuracy-weighted branching efficiency
├── fitBeta()               # Log-linear β regression over depth samples
│
├── bfsExpand()             # BFS implementation with EDBT instrumentation
├── dijkstraExpand()        # Dijkstra with EDBT instrumentation
├── astarExpand()           # A* / WA* with configurable weight
├── gbfsExpand()            # Greedy BFS
├── biDijkstraExpand()      # Bidirectional Dijkstra
├── rbfsExpand()            # RBFS with IDA*-style f-limit escalation
├── idaExpand()             # IDA* with EDBT instrumentation
├── abExpand()              # Alpha-Beta on synthetic game tree
│
├── runExperiments()        # Main experiment loop across all depths × trials
├── printBetaTable() etc.   # Console report printers
│
├── drawLineChart()         # Swing: node expansions line chart (log scale)
└── drawBEIBars()           # Swing: ABEI bar chart
```

---

## Metrics Explained

| Metric | Symbol | Description |
|---|---|---|
| Mean Nodes | N̄ | Average node expansions per trial |
| Useful Nodes | U | Path length (nodes on solution) |
| Path Accuracy | Δᵈ | Product of per-step Δˢ along solution path |
| Step Accuracy | δs | Mean Δˢ across all expansion steps |
| Pruning Ratio | Δᶜ | Mean fraction of successors discarded |
| BEI | — | Branching efficiency without accuracy weighting |
| ABEI | — | Accuracy-weighted branching efficiency |
| ACP | N̄·Δᵈ | Accuracy-Complexity Product, verified vs γ⁻ |
| q\* | — | Convergence quality bound |

---

## License

MIT — free to use, modify, and distribute.
