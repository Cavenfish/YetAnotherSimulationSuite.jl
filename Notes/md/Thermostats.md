
All MD temperatures are calculated via

$$
T_{MD} = \frac{2 E_{kin}}{3 N_f k_{B}}
$$
where $N_f$ is the number of degrees of freedom (typically $N_f = 3 N - 3$ with $N$ being the number of atoms). The minus 3 comes from the center of mass restriction. 

## Langevin

For this thermostat the Equations of Motion are revised as 

$$
m_i \ddot{x}_i =  - \nabla V - \gamma m_i \dot{x}_i + \sqrt{2 \gamma m_i k_B T} f_i
$$

where $f_i$ is a normal distribution with $\sigma =  1$ and $\mu = 0$ . Thus, the effect on the MD simulation is the following acceleration change.

$$
\Delta \ddot{x} =  -\gamma \dot{x} + \frac{1}{m}\sqrt{2 \gamma m_i k_B T} f_i
$$
The damping parameters $\gamma$ needs to be picked carefully.


## Canonical Sampling through Velocity Rescaling 

This method was taken from [this paper](https://pubs.aip.org/aip/jcp/article-abstract/126/1/014101/186581/Canonical-sampling-through-velocity-rescaling?redirectedFrom=fulltext), but note that it has a typo and so the original source code is more instructive. Here I will just outline the math that is actually implemented in the code (namely the appendix of the manuscript).


First, if the current simulation temperature ($T_{sim}$) is zero, then we do nothing (since continuing would result in division by zero). This lets the system naturally develop some motion from its potential (or numerical noise for perfectly optimized systems).

This thermostat is a velocity scaling method, so we at each integration we scale the velocity ($v$) by $\alpha$ to get the rescaled velocities ($v'$). The scaling factor is given by 

$$
\alpha = \sqrt{\frac{K'}{K}}
$$
where $K$ is the kinetic energy and $K'$ is the new desired kinetic energy for the desired temperature. To be explicit, primes denote properties based on the selected temperature for the thermostat ($T'$).

$K'$ is given by (in the same format as the code)

$$
K' = K + (1 - c_1) (\sigma (R_2 + R_1^2) \frac{1}{N_f} - K) + 2 R_1 \sqrt{K \sigma \frac{1}{N_f} (1 - c_1) c_1}
$$

where 

$$
\begin{align}
c_1 &= e^{\frac{-1}{\tau}} \\
\sigma &= \frac{1}{2} N_f k_B T' \\
R_1 &= rand(Normal(\mu=0.0, \sigma=1.0)) \\
\text{if $N_f - 1$ is odd:}\\
R_2 &= 2 * rand(Gamma(\alpha=\frac{N_f -2 }{2}, \theta=1)) + rand(Normal(\mu=0.0, \sigma=1.0))^2 \\
\text{if $N_f - 1$ is even:} \\
R_2 &= 2 * rand(Gamma(\alpha=\frac{N_f -1 }{2}, \theta = 1))
\end{align}
$$

All random selections from distributions (Normal and Gamma) are done using the [Distributions.jl](https://juliastats.org/Distributions.jl/stable/) package.