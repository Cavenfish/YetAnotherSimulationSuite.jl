Within **JMD** this potential is titled HGNN and the original manuscript is [Chen et al. 2020.](https://pubs.aip.org/aip/jcp/article/153/5/054310/1065758/Energy-transfer-between-vibrationally-excited)

This is a permutation invariant polynomial (PIP) neural network (NN) based potential. The symmetry functions (or PIP functions) used are

$$
\begin{align}
G_1 &= p_{12} + p_{34}\\
G_2 &= p_{14} + p_{32}\\
G_3 &= \sqrt{p_{12}^2 + p_{34}^2}\\
G_4 &= \sqrt{p_{14}^2 + p_{32}^2}\\
G_5 &= \sqrt{p_{12}p_{14} + p_{32}p_{34}}\\
G_6 &= p_{13}\\
G_7 &= p_{24}
\end{align}
$$
where $p_{ij}$ is given by 

$$
p_{ij} = \text{exp}(-\eta r_{ij})
$$

with $\eta$ (0.3 in their model) being a constant and $r_{ij}$ being the distance between atoms $i$ and $j$. Note, for each pair of CO atoms the indices go as follows: $O_1, C_1, O_2, C_2 = 1,2,3,4$. 

The energy of the system is given by $E  = E_\text{inter} + E_\text{intra}$, where 

$$
E_\text{intra} = b_2 + \vec{w}_2 \cdot f_1(\vec{b}_1 + \vec{w}_1 x)
$$
with $f_1 = tanh$, $b_n$ & $w_n$ are the weights and biases from the NN. The $x$ term is a mapped version of the bond length, given by

$$
x = \frac{2(r-r_a)}{(r_b-r_a)} - 1
$$
where $r$ is the bond length,  and $r_a$ & $r_b$ are constants. 

Then the inter-molecular energy is given by

$$
E_\text{inter} = b_3 + \vec{w}_3 \cdot f_3(\vec{b}_2 + W_2 f_1(\vec{b}_1 + W_1 \vec{G}))
$$


