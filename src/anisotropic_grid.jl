function anisotropic_grid(gridedData, params, dcoffset, minVal, maxVal)
    
    # initialise useful variables
    lengthX = gridedData.lenX
    lengthY = gridedData.lenY
    gridedDataDerivIterate = init_cell(lengthX, lengthY)
    realReplace = init_cell(lengthX, lengthY)
    previousGrid = init_cell(lengthX, lengthY)
    looping = 0
    currentLoop = 0
    previousMean = 0.0
    previousMedian = 0.0
    numStops = 0
    meandiffs = zeros(Float64, params["maxLoop"])

    while looping == 0
        # ******************************
        # DerivV3
        # ******************************

        # this is the first loop, therefore need the data from gridedData
        if (currentLoop == 0)
            copy_cell(gridedData, gridedDataDerivIterate)
            print("  Starting anisotropic gridding loop, loop counter: ")
        # this is NOT the first loop, therefore need the data from realReplace
        else
            copy_cell(realReplace, gridedDataDerivIterate)
            copy_cell(realReplace, previousGrid)
        end
    
        gridedData, gridedDataDerivIterate = taylor_estimation(gridedData, gridedDataDerivIterate, params, missingdata)

        # ******************************
        # ANISOTROPY CALCULATOR
        # ******************************

        # first smooth the entire grid to remove noise issues
        smoothGrid = smoothgrid(gridedData, gridedDataDerivIterate, missingdata)

        # now calculate the structure tensor and get the eigenvector for directionality information
        m1, v11, v12 = structuretensor(gridedDataDerivIterate, smoothGrid)


        # ******************************
        # TRENDING DETERMINATION
        # ******************************
        trendLOG = trendcalc(gridedDataDerivIterate, m1, params)

        # ******************************
        # REAL DATA SCALING
        # ******************************
        multiplierCells = zeros(Float64, lengthX, lengthY)
        iMultCells = zeros(Float64, lengthX, lengthY)

        # find the real data, and the multiplier for each real data cell
        for i in 1:lengthX
            for j in 1:lengthY
                if gridedData.Flag[i,j] == 1
                    tempMult = abs(gridedData.Value[i,j] / gridedDataDerivIterate.Value[i,j])
                    multiplierCells[i,j] = tempMult
                    iMultCells[i,j] = tempMult
                end
            end
        end

        for i in 1:lengthX
            for j in 1:lengthY
                realReplace.X[i,j] = gridedDataDerivIterate.X[i,j]
                realReplace.Y[i,j] = gridedDataDerivIterate.Y[i,j]

                #  if we have data in the cell, copy to realReplace
                if (gridedData.Flag[i,j] == 1)
                    realReplace.Value[i,j] = gridedData.Value[i,j]
                    realReplace.Flag[i,j] = 1
                #  if the cell is too far from data to be used, flag realReplace
                elseif (gridedDataDerivIterate.Flag[i,j] == -1)
                    realReplace.Value[i,j] = missingdata
                    realReplace.Flag[i,j] = -1
                    iMultCells[i,j] = missingdata
                # otherwise it is interpolatedΩdata,
                # and therefore we must find the real data along lines of isotropy
                else
                    failed = false
                    multFind = false
                    iMult = 0.0
                    multI = 0
                    multJ = 0
                    eigenDevi = 0.0
                    eigenSide = 1.0
                    tenD = pi * params["angleSearch"] / 180.0

                    # have we found any multiplier
                    while multFind == false
                        origXD = 0.0
                        origYD = 0.0
                        # We need to determine the direction to search. First we can just use the eigenvector result, however if that fails we must begin searching elsewhere.
                        if eigenDevi > 0.0
                            # instead let's find the direction of the eigenvector, and move the angle by theta degrees
                            # measuring from the x-axis, so range is -90 -> +90
                            origAng = atan(v12[i,j] / v11[i,j])
                            newAng = 0.0

                            # if the initial x direction is stronger, then we will first search more towards the x-axis. Also means that if it is 0 degrees, we just assume one direction.
                            if abs(v11[i,j]) > abs(v12[i,j])
                                # then we are even, and therefore have already tried the x direction first
                                if mod(eigenDevi, 2.0) == 0.0
                                    # if it's positive, then we need to add the angle to move towards the y-axis
                                    if (v11[i,j] > 0)
                                        newAng = origAng + (eigenSide * tenD)
                                    else
                                        newAng = origAng - (eigenSide * tenD)
                                    end
                                    eigenSide += 1
                                # it is odd, and it is time to check the x direction
                                else
                                    # if it's positive, then we need to subtract the angle to move towards the x-axis
                                    if v11[i,j] > 0
                                        newAng = origAng - (eigenSide * tenD)
                                    else
                                        newAng = origAng + (eigenSide * tenD)
                                    end
                                end
                                # if the initial y direction is stronger, then we will first search more towards the y-axis. Also means that if it is -90 or 90 degrees, we just assume one direction.
                            elseif abs(v11[i,j]) < abs(v12[i,j])
                                # then we are even, and therefore have already tried the y direction first
                                if mod(eigenDevi, 2.0) == 0
                                    # if it's positive, then we need to subtract the angle to move towards the x-axis
                                    if v12[i,j] > 0
                                        newAng = origAng - (eigenSide * tenD)
                                    else
                                        newAng = origAng + (eigenSide * tenD)
                                    end
                                    eigenSide += 1
                                    # it is odd, and it is time to check the y direction
                                else
                                    # if it's positive, then we need to add the angle to move towards the y-axis
                                    if v12[i,j] > 0
                                        newAng = origAng + (eigenSide * tenD)
                                    else
                                        newAng = origAng - (eigenSide * tenD)
                                    end
                                end
                                # then they are equal, and we started at 45 degrees. Arbitrarily I will say we first search more towards the x-axis. May be worthwhile to change at a later time.
                            else
                                # then we are even, and therefore have already tried the x direction first
                                if mod(eigenDevi, 2.0) == 0
                                    # if it's positive, then we need to add the angle to move towards the y-axis
                                    if v11[i,j] > 0
                                        newAng = origAng + (eigenSide * tenD)
                                    else
                                        newAng = origAng - (eigenSide * tenD)
                                    end
                                    eigenSide += 1
                                    # it is odd, and it is time to check the x direction
                                else
                                    if v11[i,j] > 0 # if it's positive, then we need to subtract the angle to move towards the x-axis
                                        newAng = origAng - (eigenSide * tenD)
                                    else
                                        newAng = origAng + (eigenSide * tenD)
                                    end
                                end
                            end

                            # now that we have a new angle, we need to give our new x and y distances
                            # then we are a stronger x direction
                            if (abs(newAng) < (pi / 4) || abs(newAng) > (3 * pi / 4))
                                # normalize such that x is 1, then multiply by 1.1
                                origXD = 1.1
                                origYD = abs(tan(newAng)) * 1.1
                                multI = 1
                            elseif (abs(newAng) == (pi / 4) || abs(newAng) == (3 * pi / 4))
                                origXD = 1.1
                                origYD = 1.1
                                multI = 1
                            else
                                # normalize such that y is 1, then multiply by 1.1
                                origXD = abs((1 / tan(newAng))) * 1.1
                                origYD = 1.1
                                multI = 1
                            end
                            if newAng < 0
                                multJ = -1
                            else
                                multJ = 1
                            end
                        else
                            # normalize such that the large of the two eigenvectors is 1, then multiply by 1.1
                            if abs(v11[i,j]) > abs(v12[i,j])
                                # begin at 1.1 cells
                                origXD = 1.1
                                origYD = (abs(v12[i,j]) / abs(v11[i,j])) * 1.1
                            else
                                origXD = (abs(v11[i,j]) / abs(v12[i,j])) * 1.1
                                # begin at 1.1 cells
                                origYD = 1.1
                            end
                            multI = sign(v11[i,j])
                            multJ = sign(v12[i,j])
                        end

                        # Now with the x and y found, we need to look in both positive and negative directions,
                        # and find the locations of the first real data cells we run into.
                        # start at the point where it would be in the next cell if looking directly along axes
                        searchStep = 1.0
                        searching = true
                        searchPos = false
                        posGood = false
                        searchNeg = false
                        negGood = false
                        iT1 = 0
                        jT1 = 0
                        iT1_2 = 0
                        jT1_2 = 0
                        iT1_3 = 0
                        jT1_3 = 0
                        iT2 = 0
                        jT2 = 0
                        iT2_2 = 0
                        jT2_2 = 0
                        iT2_3 = 0
                        jT2_3 = 0
                        foundNearbyCellPos = false
                        foundNearbyCellPos2 = false
                        foundNearbyCellNeg = false
                        foundNearbyCellNeg2 = false
                        # begin searching with current direction
                        while searching == true
                            # find the distance we are looking
                            currentStep = 1.0 + (params["searchStepSize"] * searchStep)
                            tempXD = origXD * currentStep
                            tempYD = origYD * currentStep

                            # now we need to know which direction we are looking, and add the current distance to i and j (current location)
                            # we have not found a real data cell in the positive direction yet
                            if searchPos == false
                                tempI = Int(multI * floor(tempXD)) + 1 # + 1
                                tempJ = Int(multJ * floor(tempYD)) + 1 # + 1
                                # hit an edge in x-direction
                                if (i + tempI >= lengthX || i + tempI < 1)
                                    searchPos = true
                                    posGood = false
                                # hit an edge y-direction
                                elseif (j + tempJ >= lengthY || j + tempJ < 1)
                                    searchPos = true
                                    posGood = false
                                # outside of interpolation distance
                                elseif (gridedDataDerivIterate.Flag[i + tempI,j + tempJ] == -1)
                                    searchPos = true
                                    posGood = false
                                # if true, then we have found a real data cell
                                elseif (gridedData.Flag[i + tempI,j + tempJ] == 1)
                                    # find out which side we came in from
                                    XPM = multI * tempXD - tempI + 0.5
                                    YPM = multJ * tempYD - tempJ + 0.5

                                    for k in 1:8
                                        # this cell is real
                                        if (gridedData.Flag[gridedData.close1[i + tempI,j + tempJ,k,1],gridedData.close1[i + tempI,j + tempJ,k,2]] == 1)
                                            # and not a repeat value
                                            if (gridedData.close1[i + tempI,j + tempJ,k,1] == (i + tempI) && gridedData.close1[i + tempI,j + tempJ,k,2] == (j + tempJ))
                                            else
                                                wayX = 0
                                                wayY = 0
                                                # basically just determine which direction this is. Not the best way to do it, but should work.
                                                if (k == 7 || k == 0 || k == 1)
                                                    wayX = -1
                                                elseif (k == 3 || k == 4 || k == 5)
                                                    wayX = 1
                                                end
                                                if (k == 1 || k == 2 || k == 3)
                                                    wayY = 1
                                                elseif (k == 5 || k == 6 || k == 7)
                                                    wayY = -1
                                                end
                                                if (sign(wayX) == sign(XPM) || wayX == 0)
                                                    if (sign(wayY) == sign(YPM) || wayY == 0)
                                                        # in the rare case where two are found, we need to average over all
                                                        if foundNearbyCellPos == true
                                                            iT1_3 = tempI + wayX
                                                            jT1_3 = tempJ + wayY
                                                            foundNearbyCellPos2 = true
                                                        else
                                                            iT1_2 = tempI + wayX
                                                            jT1_2 = tempJ + wayY + 1 # + 1
                                                            foundNearbyCellPos = true
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                    iT1 = tempI + 1 # + 1
                                    jT1 = tempJ + 1 # + 1
                                    searchPos = true
                                    posGood = true
                                end
                            end
                            # we have not found a real data cell in the negative direction yet
                            if searchNeg == false
                                tempI = Int(-multI * floor(tempXD)) + 1 # + 1
                                tempJ = Int(-multJ * floor(tempYD)) + 1 # + 1
                                # hit an edge in x-direction
                                if (i + tempI >= lengthX || i + tempI < 1)
                                    searchNeg = true
                                    negGood = false
                                    # hit an edge y-direction
                                elseif (j + tempJ >= lengthY || j + tempJ < 1)
                                    searchNeg = true
                                    negGood = false
                                    # outside of interpolation distance
                                elseif (gridedDataDerivIterate.Flag[i + tempI,j + tempJ] == -1)
                                    searchNeg = true
                                    negGood = false
                                    # if true, then we have found a real data cell
                                elseif (gridedData.Flag[i + tempI,j + tempJ] == 1)
                                    # find out which side we came in from
                                    XPM = 0.5 - multI * tempXD - tempI
                                    YPM = 0.5 - multJ * tempYD - tempJ

                                    for k in 1:8
                                        # this cell is real
                                        if (gridedData.Flag[gridedData.close1[i + tempI,j + tempJ,k,1],gridedData.close1[i + tempI,j + tempJ,k,2]] == 1)
                                            # and not a repeat value
                                            if (gridedData.close1[i + tempI,j + tempJ,k,1] == (i + tempI) && gridedData.close1[i + tempI,j + tempJ,k,2] == (j + tempJ))
                                            else
                                                wayX = 0
                                                wayY = 0
                                                # basically just determine which direction this is. Not the best way to do it, but should work.
                                                if (k == 7 || k == 0 || k == 1)
                                                    wayX = -1
                                                elseif (k == 3 || k == 4 || k == 5)
                                                    wayX = 1
                                                end
                                                if (k == 1 || k == 2 || k == 3)
                                                    wayY = 1
                                                elseif (k == 5 || k == 6 || k == 7)
                                                    wayY = -1
                                                end
                                                if (sign(wayX) == sign(XPM) || wayX == 0)
                                                    if (sign(wayY) == sign(YPM) || wayY == 0)
                                                        # in the rare case where two are found, we need to average over all
                                                        if foundNearbyCellNeg == true
                                                            iT2_3 = tempI + wayX
                                                            jT2_3 = tempJ + wayY
                                                            foundNearbyCellNeg2 = true
                                                        else
                                                            iT2_2 = tempI + wayX
                                                            jT2_2 = tempJ + wayY + 1 # + 1
                                                            foundNearbyCellNeg = true
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                    iT2 = tempI
                                    jT2 = tempJ
                                    searchNeg = true
                                    negGood = true
                                end
                            end
                            # have not yet found a real data cell in one or both directions yet
                            if (searchPos == false || searchNeg == false)
                                searchStep += 1
                                # then we have not found any real data cell within the user-defined interpolation distance
                                if (params["searchStepSize"] * searchStep >= (params["interpDist"] / params["cellSize"]) - 1)
                                    # exit the loop
                                    searching = false
                                end
                                # both are true, and therefore we are done searching
                            else
                                searching = false
                            end
                        end
                        # then both directions succeeeded
                        # if j + jT2 < 1
                        #     println("j ", j, " jT2 ", jT2)
                        # end
                        # if j + jT2_2 < 1
                        #     println("j ", j, " jT2_2 ", jT2_2)
                        # end
                        if posGood && negGood
                            multiplierCPos = 0.0
                            multiplierCNeg = 0.0
                            if foundNearbyCellPos
                                multiplierCPos = (multiplierCells[i + iT1,j + jT1] + multiplierCells[i + iT1_2,j + jT1_2]) / 2
                            elseif foundNearbyCellPos2
                                multiplierCPos = (multiplierCells[i + iT1,j + jT1] + multiplierCells[i + iT1_2,j + jT1_2] + multiplierCells[j + jT1_3,i + iT1_3]) / 3
                            else
                                multiplierCPos = multiplierCells[i + iT1,j + jT1]
                            end
                            if foundNearbyCellNeg
                                multiplierCNeg = (multiplierCells[i + iT2,j + jT2] + multiplierCells[i + iT2_2,j + jT2_2]) / 2
                            elseif foundNearbyCellNeg2
                                multiplierCNeg = (multiplierCells[i + iT2,j + jT2] + multiplierCells[i + iT2_2,j + jT1_2] + multiplierCells[i + iT2_3,j + jT2_3]) / 3
                            else
                                multiplierCNeg = multiplierCells[i + iT2,j + jT2]
                            end
                            distance1 = sqrt(iT1 * iT1 + jT1 * jT1)
                            distance2 = sqrt(iT2 * iT2 + jT2 * jT2)
                            totDist = distance1 + distance2
                            iMultT1 = distance2 * (multiplierCPos) / totDist
                            iMultT2 = distance1 * (multiplierCNeg) / totDist
                            iMult = iMultT1 + iMultT2
                            multFind = true
                        elseif posGood
                            if i + iT1 < 1 || j + jT1 < 1 || i + iT1_2 < 1 || j + jT1_2 < 1
                                println("TOO SMALL: i ", i, ", iT1 ", iT1, ", j ", j, ", jT1 ", jT1, ", iT1_2 ", iT1_2, ", jT1_2 ", jT1_2)
                            end
                            if i + iT1 > lengthX || j + jT1  > lengthY || i + iT1_2  > lengthX || j + jT1_2  > lengthY
                                println("TOO BIG: i ", i, ", iT1 ", iT1, ", j ", j, ", jT1 ", jT1, ", iT1_2 ", iT1_2, ", jT1_2 ", jT1_2)
                            end
                            if foundNearbyCellPos
                                iMult = (multiplierCells[i + iT1,j + jT1] + multiplierCells[i + iT1_2,j + jT1_2]) / 2
                            else
                                iMult = multiplierCells[i + iT1,j + jT1]
                            end
                            multFind = true
                        elseif negGood
                            if foundNearbyCellNeg
                                iMult = (multiplierCells[i + iT2,j + jT2] + multiplierCells[i + iT2_2,j + jT2_2]) / 2
                            else
                                iMult = multiplierCells[i + iT2,j + jT2]
                            end
                            multFind = true
                            # both directions failed. We must now start over with a new direction.
                        else
                            # this direction did not work, so search new angle
                            eigenDevi += 1
                            # then we've checked every direction
                            if ((eigenSide + 1) * params["angleSearch"] >= 180)
                                failed = true
                                multFind = true
                            end
                        end
                    end
                    if failed
                        # failed data
                        realReplace.Value[i,j] = missingdata
                        # outside of grid
                        realReplace.Flag[i,j] = -1
                        iMultCells[i,j] = missingdata
                    else
                        iMultCells[i,j] = iMult
                    end
                end
            end
        end

        # Smooth multiplier grid
        iMultS = smoothmultiplier(gridedData, gridedDataDerivIterate, iMultCells, params, missingdata)

        # Apply the smoothed normalization
        for i in 1:lengthX
            for j in 1:lengthY
                if realReplace.Flag[i,j] != 1 && realReplace.Flag[i,j] != -1
                    tempMR = gridedDataDerivIterate.Value[i,j] + (((gridedDataDerivIterate.Value[i,j] * iMultS[i,j]) - gridedDataDerivIterate.Value[i,j]) * trendLOG[i,j])
                    # do not let the scaled data become greater than the max value found in the raw data
                    if tempMR > maxVal
                        tempMR = maxVal
                    elseif tempMR < minVal + dcoffset # do not let the scaled data become less than the min value (+ DC offset) found in the raw data
                        tempMR = minVal + dcoffset
                    end
                    # replace with scaled data
                    realReplace.Value[i,j] = tempMR
                    # not real
                    realReplace.Flag[i,j] = 0
                end
            end
        end
        # ***********************************

        # Stopping criteria check************
        if params["autoStop"]
            if currentLoop != 0
                diffs = Vector{Float64}()

                for i in 1:lengthX
                    for j in 1:lengthY
                        if realReplace.Flag[i,j] != -1
                            tempVal = abs(realReplace.Value[i,j] - previousGrid.Value[i,j])
                            push!(diffs, tempVal)
                        end
                    end
                end
                diffAvg = mean(diffs)
                #  stdev seems not to be used as a stopping criterion
                # diffStd = std(diffs)

                sort!(diffs)
                diffMedi = median(diffs)

                meandiffs[currentLoop] = diffAvg

                if currentLoop != 1
                    meanDiff = diffAvg - previousMean
                    medianDiff = diffMedi - previousMedian

                    if (meanDiff >= 0 && medianDiff >= 0)
                        numStops += 1
                    else
                        previousMean = diffAvg
                        previousMedian = diffMedi
                    end
                else
                    previousMean = diffAvg
                    previousMedian = diffMedi
                end
            end
        end
        currentLoop += 1
        if (params["autoStop"] == true && numStops >= 3)
            looping = 1
        end
        # Hard stop of max iterations
        if currentLoop == params["maxLoop"]
            looping = 1
        end
        print(" ", currentLoop)
    end
    return realReplace
end
