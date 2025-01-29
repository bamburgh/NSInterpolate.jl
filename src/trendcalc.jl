# *************TRENDING DETERMINATION
function trendcalc(gridedDataDerivIterate, m1, params)
    lengthX = gridedDataDerivIterate.lenX
    lengthY = gridedDataDerivIterate.lenY

    eigenAvg = 0.0
    flatEigenVal = Vector{Float64}()

    for i in 1:lengthX
        for j in 1:lengthY
            if (gridedDataDerivIterate.Flag[i,j] != -1)
                push!(flatEigenVal, m1[i,j])
            end
        end
    end
    # sort the list
    sort!(flatEigenVal)

    # trendM larger = less change during normalization
    if params["trendM"] == 0
        eigenAvg = flatEigenVal[length(flatEigenVal) - 1]
    else
        # params["trendM"] larger = less change during normalization
        # therefore to make this work we take 100 - params["trendM"]
        tempTrendM = convert(Int64, ceil((100.0 - params["trendM"]) / 100.0))
        eigenAvg = flatEigenVal[length(flatEigenVal) * tempTrendM]
    end

    # Create "Trending" grid
    trendLOG = zeros(Float64, lengthX, lengthY)
    for i in 1:lengthX
        for j in 1:lengthY
            # If the data is outside of the grid, then don't use it
            if (gridedDataDerivIterate.Flag[i,j] == -1)
                trendLOG[i,j] = 0.0
            else
                if (m1[i,j] >= eigenAvg)
                    trendLOG[i,j] = 1.0
                else
                    # the larger m1 is, the smaller trendLOG will be, and the less change that will occur during normalization
                    trendLOG[i,j] = m1[i,j] / eigenAvg
                end
            end
        end
    end
    return trendLOG
end
