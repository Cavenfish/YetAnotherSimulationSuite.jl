
function savGol(y, windowSize, order)

  #half window size
  hws = div(windowSize, 2) 

  #build z vector
  z = collect(-hws:hws)

  #Build Vandermonde matrix
  J = zeros(order, order)

  #Fill matrix
  for i in 1:order
    @. J[:, i] = z ^ (i-1)
  end
  
  #solve for coefficients (adjoint here is the same as transpose)
  a = inv(J'J) * J' * y
  
  

end