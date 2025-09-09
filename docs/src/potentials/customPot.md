# Custom Potentials

`YASS.jl` allows you to implement custom potentials for molecular simulations. This guide demonstrates how to create a custom potential using a Lennard-Jones potential for gold (Au) as an example. To create a custom potential, you need three components:

1. `PotVars` struct - Holds potential parameters
2. Initializer function - Sets up the potential
3. Evaluation functions - Calculate energies and forces

```julia
using YASS

struct AuLJPotVars{F<:Float64} <: YASS.PotVars
    ε::F    # Well depth (eV)
    σ::F    # Equilibrium distance (Å)
    rc::F   # Cutoff radius (Å)
    rc2::F  # Squared cutoff (Å²)
end

# Initialize potential with Au parameters
function AuLJPotential(x::Union{Vector{YASS.MyAtoms}, YASS.MyCell})
    ε = 5.29e-3  # eV
    σ = 2.951    # Å (minimum at 2^(1/6)σ ≈ 3.31Å)
    rc = 8.0     # Å
    rc2 = rc^2   # Precalculate squared cutoff
    
    AuLJPotVars(ε, σ, rc, rc2)
end

# Energy-only evaluation
function AuLJEnergy(u, vars)
    E = 0.0
    potVars = vars.potVars
    
    # Extract parameters
    ε = potVars.ε
    σ = potVars.σ
    rc2 = potVars.rc2
    
    # Calculate pair interactions
    for i = 1:length(u)
        for j = i+1:length(u)
            
            # Get distance vector and magnitude squared
            r = u[j] - u[i]
            r2 = dot(r, r)
            
            # Skip if beyond cutoff
            r2 > rc2 && continue
            
            # Calculate LJ terms
            σ_r2 = (σ^2)/r2
            σ_r6 = σ_r2^3
            σ_r12 = σ_r6^2
            
            # Add pair energy
            E += 4ε*(σ_r12 - σ_r6)
        end
    end
    
    E
end

# Force-only evaluation (in-place)
function AuLJForces!(F, u, vars)
    potVars = vars.potVars
    
    # Extract parameters
    ε = potVars.ε
    σ = potVars.σ
    rc2 = potVars.rc2
    
    # Calculate forces between pairs
    for i = 1:length(u)
        for j = i+1:length(u)
            
            # Get distance vector and magnitude squared
            r = u[j] - u[i]
            r2 = dot(r, r)
            
            # Skip if beyond cutoff
            r2 > rc2 && continue
            
            # Calculate LJ terms
            σ_r2 = (σ^2)/r2
            σ_r6 = σ_r2^3
            σ_r12 = σ_r6^2
            
            # Calculate force (negative gradient)
            f = 24ε/r2 * (2σ_r12 - σ_r6) .* r
            
            # Add forces
            F[i] .-= f
            F[j] .+= f
        end
    end
end

# Combined energy and forces evaluation
function AuLJEnergyForces!(F, u, vars)
    E = 0.0
    potVars = vars.potVars
    
    # Extract parameters
    ε = potVars.ε
    σ = potVars.σ
    rc2 = potVars.rc2
    
    # Calculate interactions
    for i = 1:length(u)
        for j = i+1:length(u)
            
            # Get distance vector and magnitude squared
            r = u[j] - u[i]
            r2 = dot(r, r)
            
            # Skip if beyond cutoff
            r2 > rc2 && continue
            
            # Calculate LJ terms
            σ_r2 = (σ^2)/r2
            σ_r6 = σ_r2^3
            σ_r12 = σ_r6^2
            
            # Add pair energy
            E += 4ε*(σ_r12 - σ_r6)
            
            # Calculate and add forces
            f = 24ε/r2 * (2σ_r12 - σ_r6) .* r
            F[i] .-= f
            F[j] .+= f
        end
    end
    
    E
end

# Initializer function
function AuLJ(; constraints=nothing)
    Calculator(
        AuLJPotential;
        E = AuLJEnergy,
        F = AuLJForces!,
        EF = AuLJEnergyForces!,
        constraints=constraints
    )
end
```


# Usage example

```julia
# Create gold cluster
atoms::Vector{YASS.MyAtoms} = [
    Particle([0.0, 0.0, 0.0], zeros(3), 196.97, "Au"),
    Particle([3.0, 0.0, 0.0], zeros(3), 196.97, "Au"),
    Particle([0.0, 3.0, 0.0], zeros(3), 196.97, "Au")
]

# Create calculator
calc = AuLJ()

# Get energy
E = getPotEnergy(calc, atoms)
println("Potential energy: $E eV")

# Get forces
F = getForces(calc, atoms)
println("Forces on first atom: $(F[1]) eV/Å")

# Run MD simulation
traj = run(calc, atoms, (0.0, 10.0ps), 1.0fs, NVE())
println("Trajectory length: $(length(traj)) frames")
```

## Key Points About Custom Potentials

1. Parameter Storage
   - Use `PotVars` to store constants and parameters
   - Precalculate frequently used values (like rc²)
   - Include units in comments for clarity

2. Force Calculations
   - Use in-place operations with `F[i] .-= f` style updates

3. Optimization Tips
   - Implement distance cutoffs for efficiency
   - Precalculate squared terms where possible
   - Use inplace operations and buffers to reduce memory allocations

4. Good Practices
   - Include references for potential parameters
   - Document units clearly
   - Implement all three evaluation functions for flexibility
   - Test energy conservation in MD simulations

For more examples, check the source code of built-in potentials in `YASS.jl`.
