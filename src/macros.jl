
# Force free var memory
macro free(x)
  quote
    $x = nothing
    GC.gc()
  end
end 