

function findPeaks(arr; min=0.0, max=1e30, width=1)
  pks = Int32[]

  for i in 2:length(arr)-1

    arr[i-1] < arr[i] > arr[i+1] || continue
    arr[i] > min || continue
    arr[i] < max || continue
    arr[i-width] < arr[i] > arr[i+width] || continue

    push!(pks, i)
  end 

  pks
end

function findTurningPoints(arr)
  pks = Int32[]

  for i in 2:length(arr)-1

    if arr[i-1] < arr[i] > arr[i+1]
      push!(pks, i)
    end

    if arr[i-1] > arr[i] < arr[i+1]
      push!(pks, i)
    end
  end 

  pks
end