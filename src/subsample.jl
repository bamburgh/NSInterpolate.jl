function subsample(X, Y, Xmax, Xmin, Ymax, Ymin, params, minVal, Value, realReplace, missingdata, dcoffset)

    # *******************SUBTRACT MIN VAL
    # Subtract off the min val found earlier, as well as the DC offset
    lengthX = realReplace.lenX
    lengthY = realReplace.lenY
    if minVal < 1
        for i in 1:lengthX
            for j in 1:lengthY
                if (realReplace.Flag[i, j] != -1)
                    realReplace.Value[i, j] = realReplace.Value[i, j] - dcoffset
                    if minVal < 0.0
                        realReplace.Value[i, j] = realReplace.Value[i, j] - abs(minVal)
                    end
                end
            end
        end
    end

    lengthXF = Int(ceil(((Xmax - Xmin) / params["cellSizeF"])))
    # if the division is perfect then we need to add one
    if mod((Xmax - Xmin), params["cellSize"]) == 0.0
        lengthXF += 1
    end
    lengthYF = Int(ceil(((Ymax - Ymin) / params["cellSizeF"])))
    # if the division is perfect then we need to add one
    if mod((Ymax - Ymin), params["cellSizeF"]) == 0
        lengthYF += 1
    end

    finalData = init_cell(lengthXF, lengthYF)
    finalPos = init_position(lengthXF, lengthYF)

    # **************************SUBSAMPLE MHD
    if params["cellSize"] != params["cellSizeF"]
        println("\nSubsampling") #  Later, just give this to ContentView to display
        # Find the edges of the total grid area

        # Make sure we only use original data at first, and get those cells completed. After that we assign any data to the remaining cells.
        for k in 1:length(X)
            # ints, will round down, therefore the xPos will be in the grid as using the bottom left point as the reference.
            xPos = Int(floor((X[k] - Xmin) / params["cellSizeF"])) + 1 # + 1
            yPos = Int(floor((Y[k] - Ymin) / params["cellSizeF"])) + 1 # + 1

            # Add the x pos of that reading to the cell, ...
            if params["realGridLocations"]
                finalData.X[xPos,yPos] += X[k] - finalPos.X[xPos,yPos]
                finalData.Y[xPos,yPos] += Y[k] - finalPos.Y[xPos,yPos]
            # ..., or round the position to the center of the cell
            else
                finalPos.X[xPos,yPos] = xPos * params["cellSizeF"] + Xmin + params["cellSizeF"] / 2.0
                finalPos.Y[xPos,yPos] = yPos * params["cellSizeF"] + Ymin + params["cellSizeF"] / 2.0
            end
            # Add the value of that reading to the cell.
            finalData.Value[xPos,yPos] += Value[k]
            # Account for how many readings have been assigned to the cell.
            finalData.Flag[xPos,yPos] += 1
        end

        # Now go through all original grid cells and assign them to new grid cells that have no data currently.
        for i in 1:lengthX#-1
            for j in 1:lengthY#-1
                # ints, will round down, therefore the xPos will be in the grid as using the bottom left point as the reference.
                xPos = Int(floor((realReplace.X[i,j] - Xmin) / params["cellSizeF"])) + 1 # + 1
                yPos = Int(floor((realReplace.Y[i,j] - Ymin) / params["cellSizeF"])) + 1 # + 1

                if xPos < 1 || xPos > lengthX || yPos < 1 || yPos > lengthY
                    println()
                    println("i j ", i, " ", j, " xPos ", xPos, " yPos ", yPos, " cellSize ", params["cellSizeF"])
                    println("realReplace.X[i,j] ", realReplace.X[i,j], " Xmin ", Xmin)
                    println("realReplace.Y[i,j] ", realReplace.Y[i,j], " Ymin ", Ymin)
                end
                # if this cell isn't a real data cell and isn't bad data, then we go ahead
                if (finalData.Flag[xPos,yPos] <= 0 && realReplace.Value[i,j] != missingdata)
                    # Essentially round the position to the center of the cell
                    finalData.X[xPos,yPos] = xPos * params["cellSizeF"] + Xmin + (params["cellSizeF"] / 2)
                    finalData.Y[xPos,yPos] = yPos * params["cellSizeF"] + Ymin + (params["cellSizeF"] / 2)
                    finalData.Value[xPos,yPos] += realReplace.Value[i,j]
                    finalData.Flag[xPos,yPos] -= 1
                end
            end
        end

        println(summary(finalData))
        # Check through all grid cells, and if a cell has more than one reading in it, average the value over the number of readings.
        for i in 1:lengthXF
            for j in 1:lengthYF
                # real data cell
                if finalData.Flag[i,j] >= 1
                    if params["realGridLocations"]
                        finalData.X[i,j] = finalData.X[i,j] / finalData.Flag[i,j] + finalPos.X[i,j]
                        finalData.Y[i,j] = finalData.Y[i,j] / finalData.Flag[i,j] + finalPos.Y[i,j]
                    end
                    finalData.Value[i,j] = finalData.Value[i,j] / (finalData.Flag[i,j])
                    # the position of real data vs interpolated data no longer matters, so just set all useable data to 1
                    finalData.Flag[i,j] = 1
                # there was no data added to this cell, and therefore it is outside of the range we can use
                elseif finalData.Flag[i,j] == 0
                    finalData.Value[i,j] = missingdata
                    finalData.Flag[i,j] = -1
                # interpolated data cell
                else
                    # must multiply by -1 to remove the negative flag values
                    finalData.Value[i,j] = -finalData.Value[i,j] / (finalData.Flag[i,j])
                    # the position of real data vs interpolated data no longer matters, so just set all useable data to 1 (leaving -1 as is)
                    finalData.Flag[i,j] = 1
                end
            end
        end
    else
        finalData = init_cell(realReplace.lenX, realReplace.lenY)
        copy_cell(realReplace, finalData)
    end
    return finalData
end
