# Introduction

!!! warning "Are you sure you want to use JMD.jl?"
    This package has been developed with the idea that no one, other than myself, would use it. A consequence of that mentality is that JMD.jl is a bit bloated with niche code, and not all general MD utilities are available. 

    Consider instead using [Molly.jl](https://juliamolsim.github.io/Molly.jl/stable/) or [NQCDynamics.jl](https://nqcd.github.io/NQCDynamics.jl/stable/).


### Installation

JMD.jl is not on the general registry, so add JMD through GitHub.

```julia-repl
pkg> add https://github.com/Cavenfish/JMD
```

#### MBX add-on

If you want to use MBX within JMD then you need to compile the MBX software into a shared object binary.

#### SCME/f add-on

If you want to use SCME/f within JMD then you need to compile and install SCME/f prior to using it.