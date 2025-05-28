# Geometry Optimizations

JMD uses [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/) for optimizations. You can optimize the geometry of molecular systems just `using JMD`. 

```julia
using JMD

# Read in bdys
bdys = readSystem("./myFile.xyz")

# Run a geometry optimization
new  = opt(TIP4P(), JMD.LBFGS(), bdys)
```

The optimization algorithms from Optim.jl are not re-exported by `JMD.jl`, for this reason you must either use `JMD.LBFGS()` (as shown above) or load `Optim.jl` within your script (as shown below). Note, for this to work you need to have already installed `Optim.jl`.

```julia
using JMD
using Optim

# Read in bdys
bdys = readSystem("./myFile.xyz")

# Run a geometry optimization
new  = opt(TIP4P(), LBFGS(), bdys)
```

Configuerable options for `Optim.jl` optimizations can be passed as keyword arguments to the `opt` function. For a full list of algorithms and configuerable options please see the `Optim.jl` documentation. Here we only give a few key optional arguments that are typically used in geometry optimizations.

  - f_abstol: Absolute tolerance in changes of the objective value. Defaults to 0.0.
  - g_abstol: Absolute tolerance in the gradient, in infinity norm. Defaults to 1e-8. For gradient free methods, this will control the main convergence tolerance, which is solver specific.
  - iterations: How many iterations will run before the algorithm gives up? Defaults to 1_000.
  - show_trace: Should a trace of the optimization algorithm's state be shown on stdout? Defaults to false.

An example of passing these optional arguments to `opt` is given below.

```julia
# You can also replace the ; with , but Julia convention is ;
new  = opt(
  TIP4P(), LBFGS(), bdys; f_abstol=1e-8, g_abstol=1e-5, 
  iterations=1000000, show_trace=true 
)
```