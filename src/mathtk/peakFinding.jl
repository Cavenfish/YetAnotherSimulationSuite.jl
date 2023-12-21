

function findPeaks(arr; min=0.0, max=1e30)
  pks = Int32[]

  for i in 2:length(arr)-1
    if arr[i-1] < arr[i] > arr[i+1] && arr[i] > min && arr[i] < max
      push!(pks, i)
    end
  end 

  pks
end