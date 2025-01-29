#    # *****************ALPHA-TRIMMED MEAN
function alpha_mean(gridedData, params)
    for dummy in 1:1
        alphaData = init_cell(gridedData.lenX, gridedData.lenY)
        for i in 1:gridedData.lenX
            for j in 1:gridedData.lenY
                alphaData.X[i,j] = gridedData.X[i,j]
                alphaData.Y[i,j] = gridedData.Y[i,j]
                alphaData.Flag[i,j] = gridedData.Flag[i,j]
            end
        end

        for i in 1:gridedData.lenX
            for j in 1:gridedData.lenY
                # If the data is outside of the grid, then don't use it
                if gridedData.Flag[i,j] != -1
                    points = zeros(Float64, 9)
                    points[1] = gridedData.Value[gridedData.close1[i,j,1,1],gridedData.close1[i,j,1,2]]
                    points[2] = gridedData.Value[gridedData.close1[i,j,2,1],gridedData.close1[i,j,2,2]]
                    points[3] = gridedData.Value[gridedData.close1[i,j,3,1],gridedData.close1[i,j,3,2]]
                    points[4] = gridedData.Value[gridedData.close1[i,j,4,1],gridedData.close1[i,j,4,2]]
                    points[5] = gridedData.Value[gridedData.close1[i,j,5,1],gridedData.close1[i,j,5,2]]
                    points[6] = gridedData.Value[gridedData.close1[i,j,6,1],gridedData.close1[i,j,6,2]]
                    points[7] = gridedData.Value[gridedData.close1[i,j,7,1],gridedData.close1[i,j,7,2]]
                    points[8] = gridedData.Value[gridedData.close1[i,j,8,1],gridedData.close1[i,j,8,2]]
                    points[9] = gridedData.Value[i,j]

                    # Make sure no data that is outside of the grid is being referenced
                    for k in 1:9
                        if points[k] == missingdata
                            points[k] = gridedData.Value[i,j]
                            gridedData.close1[i,j,k,1] = i
                            gridedData.close1[i,j,k,2] = j
                        end
                    end

                    sort!(points)
                    numAlpha = 0
                    sumAlpha = 0.0
                    for k in 3:length(points) - 2
                        sumAlpha = sumAlpha + points[k]
                        numAlpha += 1
                    end
                    alphaData.Value[i,j] = sumAlpha / numAlpha
                    difference = alphaData.Value[i,j] - gridedData.Value[i,j]
                    maxV = ceil(params["interpDist"] / params["cellSize"] / 2.0)
                    distR = maxV * 2.0
                    if gridedData.Flag[i,j] == 1
                        distR = 0
                    else
                        points2 = convert(Array{Float64, 1}, gridedData.close[i,j,:])

                        # custom min finder due to -1 values
                        for k in 1:8
                            if points2[k] != -1
                                if points2[k] < distR
                                    distR = points2[k]
                                end
                            end
                        end
                    end
                    # double newDiff = difference*(1/(1+Math.Exp(-3.5*(distR-(maxV/2)))));
                    if distR > maxV
                        distR = maxV
                    end
                    newDiff = difference * (1 / (maxV / distR))
                    if distR == 0
                        newDiff = 0
                    end
                    gridedData.Value[i,j] = gridedData.Value[i,j] + newDiff
                end
            end
        end
    end

    return gridedData
end
