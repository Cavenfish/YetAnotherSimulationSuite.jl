# Geometry Optimizations

`YASS.jl` provides geometry optimization capabilities through [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/). This section explains how to optimize molecular structures and crystal cells.

## Basic Usage

The simplest way to optimize a molecular structure is:

```julia
using YASS

# Read initial structure
molecule = readSystem("water.xyz")

# Run geometry optimization
optimized = opt(TIP4Pf(), YASS.LBFGS(), molecule)

# Save optimized structure
write("optimized.xyz", optimized)
```

## Optimization Algorithms

`YASS.jl` provides access to all Optim.jl algorithms. You can either use them through `YASS`:

```julia
# Using algorithms through YASS
opt(calc, YASS.LBFGS(), molecule)      # L-BFGS algorithm
opt(calc, YASS.ConjugateGradient(), molecule)  # Conjugate gradient
opt(calc, YASS.NelderMead(), molecule)  # Derivative-free
```

Or import Optim.jl directly:

```julia
using YASS
using Optim

# Using algorithms directly from Optim
opt(calc, LBFGS(), molecule)
opt(calc, ConjugateGradient(), molecule)
```

## Configuring Optimizations

The `opt` function accepts several keyword arguments to control the optimization:

```julia
optimized = opt(
    TIP4Pf(),           # Calculator
    LBFGS(),            # Algorithm
    molecule;           # Structure
    f_abstol=1e-8,      # Function value tolerance
    g_abstol=1e-5,      # Gradient tolerance  
    iterations=100_000, # Maximum iterations
    show_trace=true     # Show progress
)
```

Common optimization parameters:
- `f_abstol`: Tolerance for changes in energy (default: 0.0)
- `g_abstol`: Tolerance for forces/gradients (default: 1e-8)
- `iterations`: Maximum optimization steps (default: 1000)
- `show_trace`: Display optimization progress (default: false)

## Periodic Systems

For periodic systems, use a `Cell` object:

```julia
# Read periodic structure
crystal = readSystem("crystal.xyz")

# Optimize atomic positions only
optimized = opt(calc, LBFGS(), crystal)

# Optimize cell parameters (lattice vectors)
optimized = optCell(calc, LBFGS(), crystal)
```

## Constraints

You can apply constraints during optimization:

```julia
# Fix specific atoms (by index)
fixed = FixedAtoms([1,2,3])
calc = Calculator(TIP4Pf(); constraints=[fixed])

# Run constrained optimization
optimized = opt(calc, LBFGS(), molecule)
```

## Convergence Monitoring

To monitor optimization progress:

```julia
optimized = opt(
    calc, LBFGS(), molecule;
    show_trace=true,
    extended_trace=true,  # Show detailed info
    trace_simplex=true    # For Nelder-Mead
)
```

## Advanced Usage

For more control over the optimization:

```julia
# Custom convergence criteria
optimized = opt(
    calc, LBFGS(), molecule;
    x_tol=1e-6,        # Position tolerance
    f_calls_limit=1000, # Max energy evaluations
    g_calls_limit=1000  # Max gradient evaluations
)

# Use different line search method
optimized = opt(
    calc, LBFGS(), molecule;
    linesearch=LineSearches.BackTracking()
)
```

For additional options and algorithms, refer to the [Optim.jl documentation](https://julianlsolvers.github.io/Optim.jl/stable/).