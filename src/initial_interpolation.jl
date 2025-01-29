"""
Now go through all cells, and assign values to the ones that have no data currently,
or determine if they are too far away from real data to use.
During this process we will also find all closest data to each cell which will be
information needed when normalizing.

The meaning of the cell flag is changed and now becomes:
+1: a cell that has observations within its borders;
 0: a cell without observations within its borders and whose value has been interpolated;
-1: a cell too far from data to be interpolated and which will remain un-filled.
"""
function initial_interpolation(gridedData, params, Xmin, Ymin, missingdata)
    
    for i in 1:gridedData.lenX
        for j in 1:gridedData.lenY
            if gridedData.Flag[i,j] != 1
                # Essentially round the position to the center of the cell
                gridedData.X[i,j] = ((i-1) * params["cellSize"]) + Xmin + (params["cellSize"] / 2)
                gridedData.Y[i,j] = ((j-1) * params["cellSize"]) + Ymin + (params["cellSize"] / 2)
            end

            # L,UL,U,UR,R,DR,D,DL
            xLoop = [-1, -1, 0, 1, 1,  1,  0, -1]
            yLoop = [0,   1, 1, 1, 0, -1, -1, -1]
            complete = false
            closeDistVal = 0.0
            closeDist = (params["interpDist"] / params["cellSize"]) + 1.0
            closeDistNum = 0

            for k in 1:8
                loopI = 1
                complete = false

                while complete == false
                    # if (j + yLoop[k] * loopI) * (i + xLoop[k] * loopI) == 0
                    #     println(j, ' ', yLoop[k], ' ', loopI)
                    #     println(i, ' ', xLoop[k], ' ', loopI)
                    # end
                    # hit an edge in x-direction
                    if (i + xLoop[k] * loopI > gridedData.lenX || i + xLoop[k] * loopI <= 0)
                        # set that cell's close1 in that direction to its own position
                        if loopI == 1
                            gridedData.close1[i,j,k,1] = i
                            gridedData.close1[i,j,k,2] = j
                        end
                        gridedData.close[i,j,k] = -1
                        complete = true
                    # hit an edge y-direction
                    elseif (j + yLoop[k] * loopI > gridedData.lenY || j + yLoop[k] * loopI <= 0)
                        if loopI == 1
                            gridedData.close1[i,j,k,1] = i
                            gridedData.close1[i,j,k,2] = j
                        end
                        gridedData.close[i,j,k] = -1
                        complete = true
                    # outside of interpolation distance
                    # elseif gridedData.Flag[j + yLoop[k] * loopI, i + xLoop[k] * loopI] == -1
                    elseif gridedData.Flag[i + xLoop[k] * loopI, j + yLoop[k] * loopI] == -1
                        if loopI == 1
                            gridedData.close1[i,j,k,1] = i
                            gridedData.close1[i,j,k,2] = j
                        end
                        gridedData.close[i,j,k] = -1
                        complete = true
                    # everything is OK at this location
                    else
                        # if this is the first loop, indicate cell's x and y
                        if loopI == 1
                            gridedData.close1[i,j,k,1] = i + xLoop[k] * loopI
                            gridedData.close1[i,j,k,2] = j + yLoop[k] * loopI
                        end
                        # if this location is real data, then indicate distance to it
                        # if gridedData.Flag[j + yLoop[k] * loopI, i + xLoop[k] * loopI] == 1
                        if gridedData.Flag[i + xLoop[k] * loopI, j + yLoop[k] * loopI] == 1
                            gridedData.close[i,j,k] = loopI
                            complete = true
                            # distance to the current loop location
                            tempLoopI = sqrt((xLoop[k] * xLoop[k] + yLoop[k] * yLoop[k]))
                            tempLoopI = (loopI) * tempLoopI
                            # if this real data point is closer, then we use it
                            if closeDist > tempLoopI
                                closeDistNum = 1
                                closeDist = tempLoopI
                                # closeDistVal = gridedData.Value[j + yLoop[k] * loopI, i + xLoop[k] * loopI]
                                closeDistVal = gridedData.Value[i + xLoop[k] * loopI, j + yLoop[k] * loopI]
                            # if this real data point is the same distance away, then we will average over the values
                            elseif closeDist == tempLoopI
                                closeDistNum += 1
                                closeDist = tempLoopI
                                # closeDistVal = closeDistVal + gridedData.Value[j + yLoop[k] * loopI, i + xLoop[k] * loopI]
                                closeDistVal = closeDistVal + gridedData.Value[i + xLoop[k] * loopI, j + yLoop[k] * loopI]
                            end
                        # not real data, therefore must keep looking
                        else
                            loopI += 1
                            # Outside of the range of interest, therefore not found
                            if loopI > params["interpDist"] / params["cellSize"]
                                gridedData.close[i,j,k] = -1
                                complete = true
                            end
                        end
                    end
                end
                # if this was a real data cell, then we need to re-assign its closest real data as itself
                if (gridedData.Flag[i,j] == 1)
                    gridedData.close[i,j,k] = 0
                end
            end

            # if this was not a real data cell then assign the values
            if gridedData.Flag[i,j] != 1
                # Now need to find closest real data and assign its value to the current cell.
                if closeDist != (params["interpDist"] / params["cellSize"]) + 1
                    gridedData.Value[i,j] = closeDistVal / closeDistNum
                    gridedData.Flag[i,j] = 0
                # no real data was close enough, therefore it is outside of the range.
                else
                    gridedData.Value[i,j] = missingdata
                    gridedData.Flag[i,j] = -1
                    for k in 1:8
                        gridedData.close1[i,j,k,1] = i
                        gridedData.close1[i,j,k,2] = j
                    end
                end
            end
        end
    end
    return gridedData
end
