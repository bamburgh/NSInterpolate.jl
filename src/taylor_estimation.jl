function taylor_estimation(gridedData, gridedDataDerivIterate, params, missingdata)
        
    # need a temporary variable that will hold the new values during a full iteration
    lengthX = gridedData.lenX
    lengthY = gridedData.lenY
    TaylorLoopVals = zeros(Float64, lengthX, lengthY)

    for  i in 1:lengthX
        for j in 1:lengthY
            tempTL = gridedDataDerivIterate.Value[i,j]
            TaylorLoopVals[i,j] = tempTL
        end
    end

    for i in 1:lengthX
        for j in 1:lengthY
            # If the data is outside of the grid, then don't use it
            if (gridedDataDerivIterate.Flag[i,j] == -1 || gridedDataDerivIterate.Value[i,j] == missingdata)
                gridedDataDerivIterate.Flag[i,j] = -1
            else
                #  8 points around center point
                points = zeros(Float64, 9)
                points[1] = gridedDataDerivIterate.Value[gridedData.close1[i,j,1,1],gridedData.close1[i,j,1,2]]
                points[2] = gridedDataDerivIterate.Value[gridedData.close1[i,j,2,1],gridedData.close1[i,j,2,2]]
                points[3] = gridedDataDerivIterate.Value[gridedData.close1[i,j,3,1],gridedData.close1[i,j,3,2]]
                points[4] = gridedDataDerivIterate.Value[gridedData.close1[i,j,4,1],gridedData.close1[i,j,4,2]]
                points[5] = gridedDataDerivIterate.Value[gridedData.close1[i,j,5,1],gridedData.close1[i,j,5,2]]
                points[6] = gridedDataDerivIterate.Value[gridedData.close1[i,j,6,1],gridedData.close1[i,j,6,2]]
                points[7] = gridedDataDerivIterate.Value[gridedData.close1[i,j,7,1],gridedData.close1[i,j,7,2]]
                points[8] = gridedDataDerivIterate.Value[gridedData.close1[i,j,8,1],gridedData.close1[i,j,8,2]]
                points[9] = gridedDataDerivIterate.Value[i,j]

                # Make sure no data that is outside of the grid is being referenced
                for k in 1:9
                    if points[k] == missingdata
                        points[k] = gridedDataDerivIterate.Value[i,j]
                        gridedData.close1[i,j,k,1] = i
                        gridedData.close1[i,j,k,2] = j
                    end
                end

                # calculate directional derivatives of `f` (Naprstek and Smith equations 2-6)
                cellarea = params["cellSize"] * params["cellSize"]
                fx = (points[5] - points[1]) / (2 * params["cellSize"])
                fy = (points[3] - points[7]) / (2 * params["cellSize"])
                fxx = (points[5] - (2 * points[9]) + points[1]) / cellarea
                fyy = (points[3] - (2 * points[9]) + points[7]) / cellarea
                fxy = (points[4] - points[2] - points[6] + points[8]) / (4 * cellarea)

                # calculate and store each successful estimate of `f(i,j)` from nearby points
                numOfSuccess = 0
                successValues = Vector{Float64}()
                for  m in -1:1
                    for  n in -1:1
                        # center point; therefore ignore
                        if (m == 0 && n == 0)
                        # make sure we don't look outside of range
                        elseif (m + i < 1 || m + i > lengthX || n + j < 1 || n + j > lengthY)
                        else
                            # There needs to be data available within the grid
                            if (gridedDataDerivIterate.Flag[i + m, j + n] != -1)
                                if (params["spatialSmooth"] == true) # if using spatial-based smoothing, then do the standard calculation
                                    # Naprstek and Smith equation 1
                                    mhdDummy = -m * fx - n * fy - 0.5 * (fxx * m * m + (2 * m * n * fxy) + (fyy * n * n))
                                    push!(successValues, gridedDataDerivIterate.Value[i + m, j + n] + mhdDummy)
                                    numOfSuccess += 1
                                else # if not, then remove the spatial aspect. This essentially turns this step into a 3x3 alpha-trimmed mean.
                                    push!(successValues, gridedDataDerivIterate.Value[i + m, j + n])
                                    numOfSuccess += 1
                                end
                            end
                        end
                    end
                end

                # calculate the mean of the estimates as the occupied value for this iteration
                if numOfSuccess > 0
                    # 2 estimates: then we need to just do an average
                    if numOfSuccess == 2
                        TaylorLoopVals[i,j] = (successValues[1] + successValues[2]) / 2
                        gridedDataDerivIterate.Flag[i,j] = 0
                    # 1 estimate: just directly assign
                    elseif numOfSuccess == 1
                        TaylorLoopVals[i,j] = successValues[1]
                        gridedDataDerivIterate.Flag[i,j] = 0
                    # more than 2: alpha trimmed mean
                    else
                        tempQuart = numOfSuccess / 4.0
                        quarterNum = Int(ceil(tempQuart))
                        maxQuarter = length(successValues) - quarterNum
                        tempNum = 0
                        tempSum = 0.0
                        # Sort according to the value, then ...
                        sort!(successValues)
                        # ..., sum the middle two quartiles ...
                        for p in quarterNum:maxQuarter
                            tempSum += successValues[p]
                            tempNum += 1
                        end
                        # ... divide by number of used estimates to get the alpha trimmed mean.
                        TaylorLoopVals[i,j] = tempSum / tempNum
                        gridedDataDerivIterate.Flag[i,j] = 0
                    end
                # this shouldn't ever happen, but if it does, then this cell is not useable
                else
                    TaylorLoopVals[i,j] = missingdata
                    gridedDataDerivIterate.Flag[i,j] = -1
                end
            end
        end
    end

    for i in 1:lengthX
        for j in 1:lengthY
            tempVal = TaylorLoopVals[i,j]
            # replace the previous iteration with the new values
            gridedDataDerivIterate.Value[i,j] = tempVal
        end
    end
    return gridedData, gridedDataDerivIterate
end
