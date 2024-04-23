
To calculate the frequencies we need the mass scaled Hessian matrix.

$$
\hat{M}_{ij}  \hat{H}(f) = \frac{1}{\sqrt{m_im_j}} \frac{\partial^2 f}{\partial q_i \partial q_j}
$$
Since we already have the analytical gradient, we can just take the Jacobian of it.

$$
\hat{H}(f) = \hat{J}(\nabla f)^T
$$

First I calculate the mass scaling matrix by taking the product of a column vector and row vector (in the code the adjoint of the same column vector). Note, since we are have 3 dimensions, the code repeats the masses 3 times per atom but I don't in the math below.

$$
\hat{M}_{ij} = \begin{bmatrix} m_1^{-0.5} \\ \vdots \\ m_n^{-0.5} \end{bmatrix} \begin{bmatrix} m_1^{-0.5} & \dots & m_n^{-0.5} \end{bmatrix} =  \begin{bmatrix} m_1^{-1} & m_1^{-0.5}m_2^{-0.5} & \dots & m_1^{-0.5}m_n^{-0.5} \\ m_2^{-0.5}m_1^{-0.5}  & \dots \\
\vdots & \ddots  \\ m_n^{-0.5}m_1^{-0.5} & \dots & & m_n^{-1} \end{bmatrix}
$$


The Hessian is calculated by using the [FiniteDifferences.jl](https://juliadiff.org/FiniteDifferences.jl/latest/) package, where I use the central finite displacement method to fill in the Jacobian matrix (which is the Hessian since its being applied to the analytical gradient). 

From there I simply element-wise multiple $\hat{M}$ with $\hat{H}$, but the matrix needs to be symmetric and numerical noise from the finite differences method prevents that from happening. To fix this, I use the `Symmetric` type cast from the [LinearAlgebra.jl](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) package to force the matrix to be symmetric. It does this by copying the top half onto the bottom half.

Lastly, I again use the LinearAlgebra.jl package to get the eigenvalues ($\lambda_i$) and eigenvectors of $\hat{M}\hat{H}$. All that is left is to take the get the frequencies (in my preferred units) from the eigenvalues. 

$$
\nu_i = \frac{ 10^{10} \hbar}{\sqrt{e \text{~amu}}} (8065.610420) \sqrt{\lambda_i} 
$$
With $e$ being the electron charge, amu being a Dalton and $\hbar$ is the reduced Planks constant. The final number is the conversion between eV and cm$^{-1}$.  