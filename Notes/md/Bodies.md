All molecular dynamics simulations revolve around the **Atom** struct.

```julia
mutable struct Atom
  r::Vector{Float64}
  v::Vector{Float64}
  m::Float64
  s::Char
end
```

Here `r` is the position of the particle, `v` is the velocity, `m` is the mass, and `s` is the symbol (*ie.* C for carbon). 