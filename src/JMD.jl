"""
Julia Molecular Dynamics

A very minimal package for me to do MD.

author: Brian C. Ferrari
"""

module JMD

  using TOML
  using JLD2
  using FFTW
  using Libdl
  using Optim
  using LsqFit
  using MiniQhull
  using DataFrames
  using Statistics
  using Clustering
  using StaticArrays
  using LinearAlgebra
  using Distributions
  using KernelDensity
  using OrdinaryDiffEq
  using SpecialFunctions
  using FiniteDifferences


  export
    #constants
    fs, ps, ns, kB,

    #io.jl
    readASExyz, readXyz, writeXyz, writeXyzTraj,

    #helpers.jl
    getFrame, getFrame!, getLastFrame, getLastFrame!,

    #bodies.jl
    # getMols, getPairs include?

    #potentials
    MvHffCO, HGNN, MBX, SPCF, TIP4P,

    #simulation.jl
    runNVE, runNVT,

    #thermostats.jl
    Berendsen, Berendsen!, Langevin, Langevin!, BDP, BDP!, BDPnT, BDPnT!,

    #post-processing.jl
    processDynamics, processDynamics!, processTmpFiles, trackVACF,
    trackEnergyDissipation, trackAllVibEnergy, trackRadialEnergy,

    #vacf.jl
    vacfInps, VDOS, getVelMas,

    #desorb.jl

    #vibrations.jl
    getHarmonicFreqs, animateMode, getModePES, getModeInteractionPES,

    #optimizations.jl
    opt,

    #structural.jl
    rdf, adf, density,

    #decayRates.jl

    #freqShifts.jl
    getINM, getMolFreq, getAllFreqs, getFvE, getFreqCoupling,

    #participationRatio.jl
    getIPR, getPR,

    #neighbors.jl
    countNearestNeighbors,

    #vibCoup.jl
    getVibCoup,

    #alphashape.jl
    alphashape,

    #savitzkyGolay.jl
    savGol,

    #peakFinding.jl
    findPeaks, findTurningPoints,

    #anneal.jl
    anneal,

    #hitAndStick.jl
    hitAndStick, HnS
  
  #end exports

  const fs   = 0.09822694788464063 # 1fs in ASE time
  const ps   = 1000 * fs
  const ns   = 1000 * ps
  const kB   = 8.617333262e-5 # eV / K

  include("./types.jl")

  include("./io.jl")
  include("./macros.jl")
  include("./helpers.jl")

  include("./lib/MBX/libmbx.jl")

  include("./md/bodies.jl")
  include("./md/potentials/MvHffCO.jl")
  include("./md/potentials/COCOff.jl")
  include("./md/potentials/HGNN.jl")
  include("./md/potentials/TIP4P.jl")
  include("./md/potentials/SPC-F.jl")
  include("./md/potentials/CH4.jl")
  include("./md/potentials/MBX.jl")
  include("./md/simulation.jl")
  include("./md/thermostats.jl")
  include("./md/post-processing.jl")
  include("./md/potentials/funcs/intra.jl")
  include("./md/potentials/funcs/inter.jl")
  include("./md/potentials/funcs/damping.jl")
  include("./md/potentials/funcs/TTMnrg.jl")

  include("./analysis/vacf.jl")
  include("./analysis/desorb.jl")
  include("./analysis/vibrations.jl")
  include("./analysis/optimizations.jl")
  include("./analysis/structral.jl")
  include("./analysis/decayRates.jl")
  include("./analysis/freqShifts.jl")
  include("./analysis/participationRatio.jl")
  include("./analysis/neighbors.jl")
  include("./analysis/vibCoup.jl")

  include("./mathtk/alphashape.jl")
  include("./mathtk/savitzkyGolay.jl")
  include("./mathtk/peakFinding.jl")

  include("./building/anneal.jl")
  include("./building/hitAndStick.jl")

  include("./QM/orcaConfig.jl")
end # module
