
To calculate the frequencies we need the Hessian matrix.

$$
\hat{H}(f) = \frac{\partial^2 f}{\partial q_i \partial q_j}
$$
Since we already have the analytical gradient, we can just take the Jacobian of it.

$$
\hat{H}(f) = \hat{J}(\nabla f)^T
$$

