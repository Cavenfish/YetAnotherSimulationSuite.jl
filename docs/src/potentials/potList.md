# Available Potentials

YASS contains a few potentials.

## Leonard Jones Potential

```julia
# Example vars for Au
d = Dict(
  "epsilon" => 0.2297,
  "sigma" => 2.95,
  "rs" => 19.0,
  "rc" => 20.0
)

calc = LJ(d)
```

## CO Potentials

### MvHff

Based on:
van Hemert, Marc C., Junko Takahashi, and Ewine F. van Dishoeck. "Molecular
dynamics study of the photodesorption of CO ice." The Journal of Physical
Chemistry A 119.24 (2015): 6354-6369.

Link:
https://pubs.acs.org/doi/full/10.1021/acs.jpca.5b02611

```julia
calc = MvHff()
```

### HGNN

Based on:
Chen, Jun, et al. "Energy transfer between vibrationally excited carbon monoxide
based on a highly accurate six-dimensional potential energy surface." The 
Journal of Chemical Physics 153.5 (2020).

Link:
https://pubs.aip.org/aip/jcp/article/153/5/054310/1065758

```julia
calc = HGNN()
```

## H$_2$O Potentials

### TIP4P/2005f

Based on:
González, M. A., & Abascal, J. L. (2011). A flexible model for water based on
TIP4P/2005. The Journal of chemical physics, 135(22).

Link:
https://pubs.aip.org/aip/jcp/article/135/22/224516/190786

```julia
calc = TIP4Pf()
```

### Simple Point Charge Models

#### SPC/F

Based on:
Toukan, Kahled, and Aneesur Rahman. "Molecular-dynamics study of atomic
motions in water." Physical Review B 31.5 (1985): 2643.

Link:
https://journals.aps.org/prb/abstract/10.1103/PhysRevB.31.2643

```julia
calc = SPC("SPC/F")
```

#### SPC/Fd

Based on:
Dang, Liem X., and B. Montgomery Pettitt. "Simple intramolecular model
potentials for water." Journal of physical chemistry 91.12 (1987): 3349-3354.

Link:
https://doi.org/10.1021/j100296a048

```julia
calc = SPC("SPC/Fd")
```

#### SPC/Fw

Based on:
Wu, Yujie, Harald L. Tepper, and Gregory A. Voth. "Flexible simple
point-charge water model with improved liquid-state properties." The Journal
of chemical physics 124.2 (2006).

Link:
https://doi.org/10.1063/1.2136877

```julia
calc = SPC("SPC/Fw")
```
