
# Not really sure how to test this
@testset "Test SavGol filter" begin
  x = collect(1:pi/100:4pi)
  y = sin.(x)

  # add noise
  noise   = length(x) |> rand
  y_noisy =  y .+ (noise .* 0.1)

  y2 = savGol(y_noisy, 63, 8)

  @test length(y) == length(y2)
end