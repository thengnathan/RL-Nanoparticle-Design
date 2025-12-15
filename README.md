# Nanoparticle Optimization via Reinforcement Learning

## Overview

This project implements an epsilon-greedy reinforcement learning algorithm to optimize the design of superparamagnetic iron oxide nanoparticles (SPIONs) for neurotransmitter detection. The optimization focuses on maximizing the delta R2 ratio, which measures the nanoparticles' sensitivity to target molecules.

**Status:** 🚧 In Progress - Additional files and features will be added

## Background

Superparamagnetic iron oxide nanoparticles can detect neurotransmitters through changes in their magnetic relaxation properties (R2 relaxation rate). The optimal nanoparticle diameter significantly affects detection sensitivity. This project uses reinforcement learning to efficiently explore the parameter space and identify optimal designs.

## Features

- **Epsilon-greedy exploration strategy** for balancing exploitation and exploration
- **Automated diameter optimization** across candidate nanoparticle sizes
- **Monte Carlo simulation integration** for R2 relaxation calculations
- **Comprehensive result tracking** including:
  - R2 values for free and aggregated states
  - Goodness-of-fit metrics (R²)
  - Delta R2 ratios (optimization objective)
- **Visualization and data persistence** for analysis

## Files

### Core Implementation

- **`epsilon_greedy_nanoparticle-1.m`** - Main reinforcement learning implementation
  - Implements epsilon-greedy action selection
  - Manages Q-value estimation and updates
  - Orchestrates the optimization loop
  - Saves results and generates visualizations

- **`Glutamate_DrWEI_Example-1.m`** - Simulation function (referenced)
  - Computes R2 relaxation values for given nanoparticle diameters
  - Returns delta R2 ratios and fit quality metrics
  - Interfaces with Monte Carlo simulation engine

## Algorithm Details

### Epsilon-Greedy Strategy

The algorithm balances exploration and exploitation using the epsilon-greedy approach:
```
With probability (1 - ε): Select action with highest Q-value (exploitation)
With probability ε: Select random action from non-optimal set (exploration)
```

### Q-Value Updates

Q-values are updated using incremental averaging:
```
Q(a) ← Q(a) + [R - Q(a)] / k(a)
```

Where:
- `Q(a)` = estimated value of action a
- `R` = observed reward (delta R2 ratio)
- `k(a)` = number of times action a has been taken

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `epsilon` | Exploration probability | 0.25 |
| `dp_input` | Candidate NP diameters (nm) | 5-60 nm in 5nm increments |
| `numTests` | Total optimization iterations | 1 (configurable) |

## Usage

### Prerequisites

- MATLAB R2018b or later
- Required functions:
  - `Glutamate_DrWEI_Example-1.m` (or equivalent simulation function)
  - Monte Carlo simulation dependencies

### Running the Optimization

1. **Set up directory structure:**
```matlab
   % The script will prompt for output directory selection
   folder = uigetdir;
```

2. **Configure parameters** (optional):
```matlab
   epsilon = 0.25;              % Adjust exploration rate
   dp_input = (1/12:1/12:1) .* 60;  % Modify diameter range
   numTests = 100;              % Set number of iterations
```

3. **Run the script:**
```matlab
   epsilon_greedy_nanoparticle
```

### Output Structure

The script creates a timestamped directory containing:
```
/output_directory/
├── YYYYMMDD_/
│   ├── initial_results.mat      # Initialization phase results
│   ├── simulation_results.mat   # Optimization iteration results
│   └── deltaR2_ratio_result.fig # Visualization of Q-values vs diameter
```

#### Output Data Format

**initial_results** (N×6 matrix):
- Column 1: NP diameter (nm)
- Column 2: R2_free
- Column 3: R² fit quality (free state)
- Column 4: R2_agg
- Column 5: R² fit quality (aggregated state)
- Column 6: delta R2 ratio

**simulation_results** (numTests×6 matrix):
- Same structure as initial_results
- One row per optimization iteration

## Methodology

### Phase 1: Initialization
- Each candidate diameter is tested once
- Initial Q-values are set to observed delta R2 ratios
- Baseline performance is established

### Phase 2: Optimization
- For each iteration:
  1. Generate random probability p ∈ (0,1)
  2. Select action based on epsilon-greedy policy
  3. Execute simulation for selected diameter
  4. Update Q-value using incremental average
  5. Record results

### Phase 3: Analysis
- Identify optimal diameter (highest Q-value)
- Generate convergence visualization
- Save all results for further analysis

## Future Development

Planned additions to this repository:

- [ ] Different nanoparticles sizing and epsilon changes

## Research Context

This work is part of ongoing research on optimizing nanoparticle designs for neurotransmitter detection applications. The epsilon-greedy approach provides a computationally efficient method for exploring the high-dimensional parameter space while maintaining exploitation of promising candidates.

## Notes

- The simulation function name may need to be updated based on your specific implementation
- Execution time depends on the complexity of the Monte Carlo simulations
- Consider adjusting `epsilon` based on exploration-exploitation trade-offs in your specific application

## Contributing

This is an active research project with the NanoAnalytic Group at California State University, Fresno under Dr. He Wei




---
