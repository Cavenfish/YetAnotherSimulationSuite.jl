"""
YetAnotherSimulationSuite.jl (YASS)

A simulation suite for molecular dynamics in Julia. 
"""
module YetAnotherSimulationSuite

  # These are part of the Julia Standard Library
  using TOML
  using Libdl
  using Statistics
  using Serialization
  using LinearAlgebra

  # These are fully external dependencies
  using JLD2
  using FFTW
  using Optim
  using PyCall
  using Chemfiles
  using Distances
  using DataFrames
  using Clustering
  using StaticArrays
  using Distributions
  using KernelDensity
  using OrdinaryDiffEq

  import FiniteDiff: JacobianCache, finite_difference_jacobian!

  using Reexport
  @reexport using Unitful

  # Wrapper for ASE calculator with SCME/f potential
  function __init__()
    @pyinclude(joinpath(@__DIR__, "lib/SCMEf/libscmef.py"))
  end

  export

    #types.jl
    MyAtoms, MyCell, MyCalc, PotVars, ThermoVars, MyThermostat,
    MyConstraint, MyTraj,

    #io.jl
    readSystem,

    #helpers.jl
    CoM, vCoM, zeroVCoM!, getForces, reducedMass,

    #energetics.jl
    vibExcite!, transExcite!, getPotEnergy, getRotEnergy, getTransEnergy,
    getVibEnergy,

    #molsAndPairs.jl
    getMols,

    #bodies.jl
    Particle, swapIso!, pickRandomMol, centerBdys!,
    swapAtoms!, center,

    #cells.jl
    makeCell, makeBdys, getScaledPos, getPos, wrap!, makeSuperCell,
    makeSuperCell!, center!, getVolume, reorder!,

    #constraints.jl
    fixAtoms,

    #trajectories.jl
    getVelMas, getLastFrame!,

    #intra.jl
    morse, morse!, harmonicBond, harmonicBond!, harmonicBondAngle,
    harmonicBondAngle!,

    #inter.jl
    vdw, vdw!, buckingham, buckingham!, coulomb, coulomb!,

    #cutoff.jl
    switchLR!, switchSR!, switchAP!,

    #potentials
    MvHff, HGNN, MBX, SPCF, TIP4Pf, SCMEf,

    #thermostats
    Berendsen, Langevin, CVR,

    #ensembles
    NVE, NVT,

    #post-processing.jl
    processDynamics, processTmpFiles, 
    
    #vacf.jl
    vacfInps, VDOS, getDiffusionCoefficient, Hann, HannM, Welch,

    #vibrations.jl
    getHarmonicFreqs, animateMode, getModePES, getModeInteractionPES,

    #optimizations.jl
    opt, hiddenOpt, optCell,

    #distributions.jl
    rdf, density,

    #geometric_operations.jl
    translate!, rotate!, randRotate!, 

    #stress.jl
    getNumericalStress, getNumericalStressOrthogonal,

    #participationRatio.jl
    getIPR, getPR,

    #neighbors.jl
    countNearestNeighbors,

    #savitzkyGolay.jl
    savGol,

    #peakFinding.jl
    findPeaks, findTurningPoints,

    #anneal.jl
    anneal,

    #hitAndStick.jl
    hitAndStick, HnS,

    #phonopy.jl
    phonopy_addForces, phonopy_getDisplacements, phonopy_getPhonons,
    phonopy_getDisplacementsDataset, reorderPhonopyForces!,

    #libscmef.jl
    scmef_getDipole,

    #libmbx.jl
    mbx_getDipole
  
  #end exports

  include("./types.jl")

  include("./io.jl")
  include("./macros.jl")
  include("./helpers.jl")

  include("./lib/MBX/libmbx.jl")
  include("./lib/SCMEf/libscmef.jl")
  include("./lib/Phonopy/phonopy.jl")

  include("./md/cells.jl")
  include("./md/bodies.jl")
  include("./md/simulation.jl")
  include("./md/calculators.jl")
  include("./md/constraints.jl")
  include("./md/thermostats.jl")
  include("./md/trajectories.jl")
  include("./md/post-processing.jl")

  include("./md/thermostats/CVR.jl")
  include("./md/thermostats/Langevin.jl")
  include("./md/thermostats/Berendsen.jl")

  include("./md/potentials/MvHff.jl")
  include("./md/potentials/HGNN.jl")
  include("./md/potentials/TIP4P.jl")
  include("./md/potentials/SPC-F.jl")
  include("./md/potentials/MBX.jl")
  include("./md/potentials/SCMEf.jl")

  include("./md/potentials/funcs/PBC.jl")
  include("./md/potentials/funcs/intra.jl")
  include("./md/potentials/funcs/inter.jl")
  include("./md/potentials/funcs/damping.jl")
  include("./md/potentials/funcs/cutoff.jl")
  include("./md/potentials/funcs/auxiliary.jl")

  include("./analysis/vacf.jl")
  include("./analysis/vibrations.jl")
  include("./analysis/optimizations.jl")
  include("./analysis/participationRatio.jl")
  include("./analysis/energetics.jl")

  include("./structural/distributions.jl")
  include("./structural/molsAndPairs.jl")
  include("./structural/neighbors.jl")

  include("./mathtk/savitzkyGolay.jl")
  include("./mathtk/peakFinding.jl")
  include("./mathtk/stress.jl")
  include("./mathtk/basicMath.jl")
  include("./mathtk/geometric_operations.jl")

  include("./building/anneal.jl")
  include("./building/hitAndStick.jl")
end # module
