# Smooth multiplier grid
function smoothmultiplier(gridedData, gridedDataDerivIterate, iMultCells, params, missingdata)
    lengthX = gridedData.lenX
    lengthY = gridedData.lenY
    iMultS = zeros(Float64, lengthX, lengthY)
    # if no smoothing is wanted, then just pass along the multiplier grid
    if params["multiSmooth"] == 0
        # don't do anything
        iMultS = iMultCells
        # otherwise we will apply a 3x3 smoother to the grid, with a sliding effect
    else
        for i in 1:lengthX
            for j in 1:lengthY
                # If the data is outside of the grid, then don't use it
                if (gridedDataDerivIterate.Flag[i,j] == -1)
                    iMultS[i,j] = 0.0
                else                    
                    points = Vector{Float64}()
                    push!(points, iMultCells[gridedData.close1[i,j,1,1],gridedData.close1[i,j,1,2]])
                    push!(points, iMultCells[gridedData.close1[i,j,2,1],gridedData.close1[i,j,2,2]])
                    push!(points, iMultCells[gridedData.close1[i,j,3,1],gridedData.close1[i,j,3,2]])
                    push!(points, iMultCells[gridedData.close1[i,j,4,1],gridedData.close1[i,j,4,2]])
                    push!(points, iMultCells[gridedData.close1[i,j,5,1],gridedData.close1[i,j,5,2]])
                    push!(points, iMultCells[gridedData.close1[i,j,6,1],gridedData.close1[i,j,6,2]])
                    push!(points, iMultCells[gridedData.close1[i,j,7,1],gridedData.close1[i,j,7,2]])
                    push!(points, iMultCells[gridedData.close1[i,j,8,1],gridedData.close1[i,j,8,2]])

                    # Make sure no data that is outside of the grid is being referenced
                    for k in 1:8
                        if points[k] == missingdata
                            points[k] = iMultCells[i,j]
                        end
                    end

                    # now to enact the sliding scale, add 100 minus the smoother amount number of values of the original mutlitplier
                    # e.g. 1 will add 100, 100 will add only 1
                    smoother = 101 - Int(params["multiSmooth"])
                    for idummy in 1:smoother
                        push!(points, iMultCells[i,j])
                    end
                    iMultS[i,j] = mean(points)
                end
            end
        end
    end
    return iMultS
end
