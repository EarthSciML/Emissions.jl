using Emissions
using Test


@testset "findlayer tests" begin
    @test findLayer([0.0, 10.0, 20.0, 30.0, 40.0], 0.0) == 1
    @test findLayer([0.0, 10.0, 20.0, 30.0, 40.0], 9.0) == 1
    @test findLayer([0.0, 10.0, 20.0, 30.0, 40.0], 15.0) == 2
    @test findLayer([0.0, 10.0, 20.0, 30.0, 40.0], 35.0) == 4
    @test_throws ErrAboveModelTop findLayer([0.0, 10.0, 20.0, 30.0, 40.0], 45.0)
end

# ╔═╡ 4f9fd573-9a81-43e3-b9f9-a810ff84ffc4
@testset "calcDeltaH tests" begin
	temperature = [50, 10, 15, 15]
	windSpeed = [10, 12.5, 15, 14]
	sClass = [0.25, 0.75, 0.4, 0.6]
	s1 = [0.2, 0.5, 1.0, NaN]
	stackHeight = 100.0
	stackTemp = 80.0
	stackVel = 20.0
	stackDiam = 10.0
	@test calcDeltaH(1, temperature, windSpeed, sClass, s1, stackHeight, stackTemp, stackVel, stackDiam) - 26.3901 < 0.0001
	@test calcDeltaH(2, temperature, windSpeed, sClass, s1,stackHeight, stackTemp, stackVel, stackDiam) - 255.8218 < 0.0001
	@test calcDeltaH(3, temperature, windSpeed, sClass, s1, stackHeight, stackTemp, stackVel, stackDiam) - 168.4916 < 0.0001
    try
        calcDeltaH(4, temperature, windSpeed, sClass, s1, stackHeight, stackTemp, 	stackVel, stackDiam) - 168.4916 < 0.0001
    catch err
		@test err isa Exception
    	@test sprint(showerror, err) == "plume height == NaN deltaH: NaN, stackDiam: 10.0, stackVel: 20.0, windSpd: 14.0, stackTemp: 80.0, airTemp: 15, stackHeight: 100.0"
    end
end

# ╔═╡ 44cea95b-10c7-4a1d-9086-8813cfa3670e
@testset "CalcDeltaHPrecomputed tests" begin
	temperature = [80.0, 80.0, 100.0, 20.0, 20.0, 20.0, 20.0]
	windSpeed = [10.0, 10.0, 40.0, 40.0, 40.0, 40.0, 40.0]
	sClass = [1.0, 1.0, 1.0, 1.0, 1.0, 0.25, 0.25]
	s1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
	windSpeedMinusOnePointFour = [NaN, 10.0, 10.0, 1.0, 1.0, 1.0, 1.0]
	windSpeedMinusThird = [1.0, 1.0, 1.0, NaN, 1.0, 1.0, 1.0]
	windSpeedInverse = [1.0, 1.0, 1.0, 1.0, 1.0, NaN, 1.0]
	stackHeight = 100.0
	stackTemp = 100.0
	stackVel = 20.0
	stackDiam = 10.0

	@test_throws ErrorException("plumerise: momentum-dominated deltaH is NaN. stackDiam: 10.0, stackVel: 20.0, windSpeedMinusOnePointFour: NaN") calcDeltaHPrecomputed(1, temperature, windSpeed, sClass, s1, stackHeight, stackTemp, stackVel, stackDiam, windSpeedMinusOnePointFour, windSpeedMinusThird, windSpeedInverse)
	 
	@test calcDeltaHPrecomputed(2, temperature, windSpeed, sClass, s1, stackHeight, stackTemp, stackVel, stackDiam, windSpeedMinusOnePointFour, windSpeedMinusThird, windSpeedInverse) - 6628.9080 < 0.0001
	
	@test calcDeltaHPrecomputed(3, temperature, windSpeed, sClass, s1, stackHeight, stackTemp, stackVel, stackDiam, windSpeedMinusOnePointFour, windSpeedMinusThird, windSpeedInverse) - 0 < 0.0001


	@test_throws ErrorException("plumerise: stable bouyancy-dominated deltaH is NaN. F: 6537.766666666667, s1: 1.0, windSpeedMinusThird: NaN") calcDeltaHPrecomputed(4, temperature,windSpeed, sClass, s1, stackHeight, stackTemp, stackVel, stackDiam, windSpeedMinusOnePointFour, windSpeedMinusThird, windSpeedInverse)
	
	@test calcDeltaHPrecomputed(5, temperature, windSpeed, sClass, s1, stackHeight, stackTemp, stackVel, stackDiam, windSpeedMinusOnePointFour, windSpeedMinusThird, windSpeedInverse) - 542.2602 < 0.0001

	@test_throws ErrorException("plumerise: unstable bouyancy-dominated deltaH is NaN. F: 6537.766666666667, stackHeight: 100.0, windSpeedInverse: NaN") calcDeltaHPrecomputed(6, temperature, windSpeed, sClass, s1, stackHeight, stackTemp, stackVel, stackDiam, windSpeedMinusOnePointFour, windSpeedMinusThird, windSpeedInverse)
		
	@test calcDeltaHPrecomputed(7, temperature, windSpeed, sClass, s1, stackHeight, stackTemp, stackVel, stackDiam, windSpeedMinusOnePointFour, windSpeedMinusThird, windSpeedInverse) - 2981.0884 < 0.0001

end

# ╔═╡ 62de0080-44b2-4266-a15b-4469062639ee
@testset "ASME tests" begin
	temperature = [50.0, 10.0, 15.0, 15.0]
	windSpeed = [10.0, 12.5, 15.0, 14.0]
	sClass = [0.25, 0.75, 0.4, 0.6]
	s1 = [0.2, 0.5, 1.0, NaN]
	stackTemp = 100.0
	stackVel = 20.0
	stackDiam = 10.0
	layerHeights = [0.0, 10.0, 20.0, 30.0, 40]

	@test calcDeltaH(1, temperature, windSpeed, sClass, s1, 0.0, stackTemp, stackVel, stackDiam) == 0

	@test ASME(0.0, stackDiam, stackTemp, stackVel, layerHeights, temperature, windSpeed, sClass, s1) == (1, 0)

	@test_throws ErrAboveModelTop ASME(10.0, stackDiam, stackTemp, stackVel, layerHeights, temperature, windSpeed, sClass, s1)

	@test_throws ErrAboveModelTop ASME(20.0, stackDiam, stackTemp, stackVel, layerHeights, temperature, windSpeed, sClass, s1)

	@test_throws ErrAboveModelTop ASME(30.0, stackDiam, stackTemp, stackVel, layerHeights, temperature, windSpeed, sClass, s1)

	@test_throws ErrorException(string("plume height == NaN ",
			  "deltaH: NaN, stackDiam: 10.0, ",
			  "stackVel: 20.0, windSpd: 14.0, stackTemp: 100.0, ",
			  "airTemp: 15.0, stackHeight: 40.0")) ASME(40.0, stackDiam, stackTemp, stackVel, layerHeights, temperature, windSpeed, sClass, s1)

	@test_throws ErrAboveModelTop ASME(50.0, stackDiam, stackTemp, stackVel, layerHeights, temperature, windSpeed, sClass, s1)

	@test_throws ErrAboveModelTop ASME(60.0, stackDiam, stackTemp, stackVel, layerHeights, temperature, windSpeed, sClass, s1)
	
end

# ╔═╡ 41bd0f8b-d3c7-4033-bbce-b5ef61604220
@testset "ASMEPrecomputed tests" begin
	temperature = [80.0, 80.0, 100.0, 20.0, 20.0, 20.0, 20.0]
	windSpeed = [10.0, 10.0, 40.0, 40.0, 40.0, 40.0, 40.0]
	sClass = [1.0, 1.0, 1.0, 1.0, 1.0, 0.25, 0.25]
	s1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
	windSpeedMinusOnePointFour = [NaN, 10.0, 10.0, 1.0, 1.0, 1.0, 1.0]
	windSpeedMinusThird = [1.0, 1.0, 1.0, NaN, 1.0, 1.0, 1.0]
	windSpeedInverse = [1.0, 1.0, 1.0, 1.0, 1.0, NaN, 1.0]
	stackTemp = 100.0
	stackVel = 20.0
	stackDiam = 10.0
	layerHeights = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0]

	@test_throws ErrorException("plumerise: momentum-dominated deltaH is NaN. stackDiam: 10.0, stackVel: 20.0, windSpeedMinusOnePointFour: NaN") ASMEPrecomputed(0.0, stackDiam, stackTemp, stackVel, layerHeights, temperature, windSpeed, sClass, s1, windSpeedMinusOnePointFour, windSpeedMinusThird, windSpeedInverse) 

	@test_throws ErrorException("plumerise: momentum-dominated deltaH is NaN. stackDiam: 10.0, stackVel: 20.0, windSpeedMinusOnePointFour: NaN") ASMEPrecomputed(10.0, stackDiam, stackTemp, stackVel, layerHeights, temperature, windSpeed, sClass, s1, windSpeedMinusOnePointFour, windSpeedMinusThird, windSpeedInverse) 

	@test_throws ErrAboveModelTop ASMEPrecomputed(20.0, stackDiam, stackTemp, stackVel, layerHeights, temperature, windSpeed, sClass, s1, windSpeedMinusOnePointFour, windSpeedMinusThird, windSpeedInverse) 	

	@test ASMEPrecomputed(30.0, stackDiam, stackTemp, stackVel, layerHeights, temperature, windSpeed, sClass, s1, windSpeedMinusOnePointFour, windSpeedMinusThird, windSpeedInverse) == (3, 30.0) 

	@test_throws ErrorException("plumerise: stable bouyancy-dominated deltaH is NaN. F: 6537.766666666667, s1: 1.0, windSpeedMinusThird: NaN") ASMEPrecomputed(40.0, stackDiam, stackTemp, stackVel, layerHeights, temperature, windSpeed, sClass, s1, windSpeedMinusOnePointFour, windSpeedMinusThird, windSpeedInverse)

	@test_throws ErrAboveModelTop ASMEPrecomputed(50.0, stackDiam, stackTemp, stackVel, layerHeights, temperature, windSpeed, sClass, s1, windSpeedMinusOnePointFour, windSpeedMinusThird, windSpeedInverse)

	@test_throws ErrorException("plumerise: unstable bouyancy-dominated deltaH is NaN. F: 6537.766666666667, stackHeight: 60.0, windSpeedInverse: NaN") ASMEPrecomputed(60.0, stackDiam, stackTemp, stackVel, layerHeights, temperature, windSpeed, sClass, s1, windSpeedMinusOnePointFour, windSpeedMinusThird, windSpeedInverse)

	@test_throws ErrAboveModelTop ASMEPrecomputed(70.0, stackDiam, stackTemp, stackVel, layerHeights, temperature, windSpeed, sClass, s1, windSpeedMinusOnePointFour, windSpeedMinusThird, windSpeedInverse)
	
end
