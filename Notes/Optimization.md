
The current geometry optimization process is a bit convoluted but it works for some algorithms. First I flatten the positions matrix into a vector, which once received by the potential calculator is "reinflated" (I make a matrix again).  Pairs and molecules are determined from the initial input and remain static during the optimization. 


Algorithms that work:
 1. LBFGS
 2. CG


TODO:
1. Make `HGNNpot` function only calculate forces if `G != nothing` 
2. Derive analytical Hessian and include it in `HGNNpot`
3. Improve performance (check allocations and memory)
4. 
