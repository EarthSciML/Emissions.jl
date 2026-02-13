module Emissions
export g,
    ErrAboveModelTop, findLayer, calcDeltaH, ASME, calcDeltaHPrecomputed, ASMEPrecomputed

const g = 9.80665

# ErrAboveModelTop is an error that is returned when the plume is above the top model layer.
ErrAboveModelTop = ErrorException("plume rise > top of grid")

function findLayer(layerHeights::Vector{Float64}, stackHeight::Float64)
    stackLayer = searchsortedfirst(layerHeights, stackHeight)
    if stackLayer > length(layerHeights)
        throw(ErrAboveModelTop)
    end
    if stackLayer != 1
        stackLayer -= 1
    end
    return stackLayer
end

# calcDeltaH calculates plume rise (ASME, 1973).
function calcDeltaH(
        stackLayer::Int,
        temperature,
        windSpeed,
        sClass,
        s1::Vector{Float64},
        stackHeight,
        stackTemp,
        stackVel,
        stackDiam::Float64
    )
    deltaH = 0.0 # Plume rise, (m).

    airTemp = temperature[stackLayer]
    windSpd = windSpeed[stackLayer]

    if (stackTemp - airTemp) < 50.0 && stackVel > windSpd && stackVel > 10.0

        # Plume is dominated by momentum forces
        deltaH = stackDiam * (stackVel^1.4) / (windSpd^1.4)

    else # Plume is dominated by buoyancy forces
        # Bouyancy flux, m4/s3
        F = g * (stackTemp - airTemp) / stackTemp * stackVel * ((stackDiam / 2)^2)

        if sClass[stackLayer] > 0.5  # stable conditions
            deltaH = 29.0 * ((F / s1[stackLayer])^0.333333333) / (windSpd^0.333333333)

        else # unstable conditions
            deltaH = 7.4 * ((F * (stackHeight^2.0))^0.333333333) / windSpd
        end
    end
    if isnan(deltaH)
        throw(
            ErrorException(
                string(
                    "plume height == NaN ",
                    "deltaH: $(deltaH), stackDiam: $(stackDiam), ",
                    "stackVel: $(stackVel), windSpd: $(windSpd), stackTemp: $(stackTemp), ",
                    "airTemp: $(airTemp), stackHeight: $(stackHeight)"
                ),
            ),
        )
        return deltaH
    end
    return deltaH
end

# ASME takes emissions stack height(m), diameter (m), temperature (K),
# and exit velocity (m/s) and calculates the k index of the equivalent
# emissions height after accounting for plume rise.
# Additional required inputs are model layer heights (staggered grid; layerHeights [m]),
# temperature at each layer [K] (unstaggered grid),
# wind speed at each layer [m/s] (unstaggered grid),
# stability class (sClass [0 or 1], unstaggered grid),
# and stability parameter (s1 [unknown units], unstaggered grid).
# Uses the plume rise calculation: ASME (1973), as described in Sienfeld and Pandis,
# ``Atmospheric Chemistry and Physics - From Air Pollution to Climate Change
function ASME(
        stackHeight,
        stackDiam,
        stackTemp,
        stackVel::Float64,
        layerHeights,
        temperature,
        windSpeed,
        sClass,
        s1::Vector{Float64}
    )
    stackLayer = findLayer(layerHeights, stackHeight)
    deltaH = calcDeltaH(
        stackLayer,
        temperature,
        windSpeed,
        sClass,
        s1,
        stackHeight,
        stackTemp,
        stackVel,
        stackDiam
    )
    print(deltaH)
    plumeHeight = stackHeight + deltaH
    plumeLayer = findLayer(layerHeights, plumeHeight)
    return plumeLayer, plumeHeight
end

# calcDeltaHPrecomputed calculates plume rise, the same as calcDeltaH,
# (ASME, 1973), except that it uses precomputed meteorological parameters.
function calcDeltaHPrecomputed(
        stackLayer::Int,
        temperature,
        windSpeed,
        sClass,
        s1::Vector{Float64},
        stackHeight,
        stackTemp,
        stackVel,
        stackDiam::Float64,
        windSpeedMinusOnePointFour,
        windSpeedMinusThird,
        windSpeedInverse::Vector{Float64}
    )
    deltaH = 0.0 # Plume rise, (m).

    airTemp = temperature[stackLayer]
    windSpd = windSpeed[stackLayer]

    if (stackTemp - airTemp) < 50.0 && stackVel > windSpd && stackVel > 10.0

        # Plume is dominated by momentum forces
        deltaH = stackDiam * (stackVel^1.4) * windSpeedMinusOnePointFour[stackLayer]

        if isnan(deltaH)
            throw(
                ErrorException(
                    string(
                        "plumerise: momentum-dominated deltaH is NaN. ",
                        "stackDiam: $(stackDiam), stackVel: $(stackVel), windSpeedMinusOnePointFour: $(windSpeedMinusOnePointFour[stackLayer])"
                    ),
                ),
            )
            return deltaH
        end

    else  # Plume is dominated by buoyancy forces
        tempDiff = 0.0
        if stackTemp - airTemp == 0
            tempDiff = 0
        else
            tempDiff = 2 * (stackTemp - airTemp) / (stackTemp + airTemp)
        end

        # Bouyancy flux, m4/s3
        F = g * tempDiff * stackVel * ((stackDiam / 2)^2)

        if sClass[stackLayer] > 0.5 && s1[stackLayer] != 0 && F > 0 # stable conditions

            # Ideally, we would also use the inverse of S1,
            # but S1 is zero sometimes so that doesn't work well.
            deltaH = 29.0 * ((F / s1[stackLayer])^0.333333333) *
                windSpeedMinusThird[stackLayer]

            if isnan(deltaH)
                throw(
                    ErrorException(
                        string(
                            "plumerise: stable bouyancy-dominated deltaH is NaN. ",
                            "F: $(F), s1: $(s1[stackLayer]), windSpeedMinusThird: $(windSpeedMinusThird[stackLayer])"
                        ),
                    ),
                )
                return deltaH
            end

        elseif F > 0.0  # unstable conditions
            deltaH = 7.4 * ((F * (stackHeight^2.0))^0.333333333) *
                windSpeedInverse[stackLayer]

            if isnan(deltaH)
                throw(
                    ErrorException(
                        string(
                            "plumerise: unstable bouyancy-dominated deltaH is NaN. ",
                            "F: $(F), stackHeight: $(stackHeight), windSpeedInverse: $(windSpeedInverse[stackLayer])"
                        ),
                    ),
                )
                return deltaH
            end
        else
            # If F < 0, the unstable algorithm above will give an imaginary plume rise.
            deltaH = 0
        end
    end
    return deltaH
end

# ASMEPrecomputed is the same as ASME except it takes
# precomputed (averaged) meteorological parameters:
# the inverse of the stability parameter (s1Inverse [1/unknown units],
# unstaggered grid), windSpeedMinusOnePointFour [(m/s)^(-1.4)] (unstaggered grid),
# windSpeedMinusThird [(m/s)^(-1/3)] (unstaggered grid),
# and windSpeedInverse [(m/s)^(-1)] (unstaggered grid),
# Uses the plume rise calculation: ASME (1973), as described in Sienfeld and Pandis,
# ``Atmospheric Chemistry and Physics - From Air Pollution to Climate Change
function ASMEPrecomputed(
        stackHeight,
        stackDiam,
        stackTemp,
        stackVel::Float64,
        layerHeights,
        temperature,
        windSpeed,
        sClass,
        s1,
        windSpeedMinusOnePointFour,
        windSpeedMinusThird,
        windSpeedInverse::Vector{Float64}
    )
    stackLayer = findLayer(layerHeights, stackHeight)
    deltaH = calcDeltaHPrecomputed(
        stackLayer,
        temperature,
        windSpeed,
        sClass,
        s1,
        stackHeight,
        stackTemp,
        stackVel,
        stackDiam,
        windSpeedMinusOnePointFour,
        windSpeedMinusThird,
        windSpeedInverse
    )

    plumeHeight = stackHeight + deltaH
    plumeLayer = findLayer(layerHeights, plumeHeight)
    return plumeLayer, plumeHeight
end

end
