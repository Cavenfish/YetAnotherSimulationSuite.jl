"""
    findPeaks(arr; min=0.0, max=1e30, width=1)

Find indices of local maxima (peaks) in an array.

# Arguments
- `arr`: Array of numeric values to search for peaks.
- `min`: Minimum value threshold for a peak (default: 0.0).
- `max`: Maximum value threshold for a peak (default: 1e30).
- `width`: Minimum width for a peak (default: 1).

# Returns
- Vector of indices where peaks are found.
"""
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

"""
    findTurningPoints(arr)

Find indices of local maxima and minima (turning points) in an array.

# Arguments
- `arr`: Array of numeric values.

# Returns
- Vector of indices where turning points are found.
"""
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