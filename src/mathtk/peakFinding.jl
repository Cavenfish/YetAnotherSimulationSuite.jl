

function findPeaks(arr; min=0.0)
  pks = Int32[]

  for i in 2:length(arr)
    if arr[i-1] < arr[i] > arr[i+1] && arr[i] > min
      push!(pks, i)
    end
  end 

  pks
end