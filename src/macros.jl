
# Force free var memory
macro free(x)
  e = quote
    $x = nothing
    GC.gc()
  end
  esc(e)
end 