
function checkOnOff(fn)
  S  = zeros(2)
  rs = 2.0
  rc = 3.0

  fn(S, 1.0, rs, rc)
  @test S[1] == 1.0
  @test S[2] == 0.0

  fn(S, 4.0, rs, rc)
  @test S[1] == 0.0
  @test S[2] == 0.0
end

@testset "Check On/Off" begin

  checkOnOff(YASS.switchSR!)
  checkOnOff(YASS.switchLR!)
  checkOnOff(YASS.switchAP!)

end
