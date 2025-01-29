# first smooth the entire grid to remove noise issues
function smoothgrid(gridedData, gridedDataDerivIterate, missingdata)
    lengthX = gridedData.lenX
    lengthY = gridedData.lenY
    smoothGrid = zeros(Float64, lengthX, lengthY)

    for i in 1:lengthX
        for j in 1:lengthY
            # If the data is outside of the grid, then don't use it
            if (gridedDataDerivIterate.Flag[i,j] != -1)
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
                    if (points[k] == missingdata)
                        points[k] = gridedDataDerivIterate.Value[i,j]
                        gridedData.close1[i,j,k,1] = i
                        gridedData.close1[i,j,k,2] = j
                    end
                end

                sumAlpha = 0.0
                for k in 1:9
                    sumAlpha = sumAlpha + points[k]
                end
                smoothGrid[i,j] = sumAlpha / 9.0
            end
        end
    end
    return smoothGrid
end
