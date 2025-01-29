function offset_positive(gridedData)
    # Find minimum real value and add to all data to ensure RealReplaceV2 works properly
    minVal = Inf #10000000000.0
    maxVal = -Inf #-10000000000.0

    for i in 1:gridedData.lenX
        for j in 1:gridedData.lenY
            if (gridedData.Flag[i,j] == 1)
                if (gridedData.Value[i,j] < minVal)
                    minVal = gridedData.Value[i,j]
                end
                if (gridedData.Value[i,j] > maxVal)
                    maxVal = gridedData.Value[i,j]
                end
            end
        end
    end

    # Need to add some positive DC offset to ensure no division by 0 (if 0s exist in this dataset)
    dcoffset = 0.0
    if minVal < 1
        dcoffset = 100.0
        # add this to the dataset
        for i in 1:gridedData.lenX
            for j in 1:gridedData.lenY
                if gridedData.Flag[i,j] != -1
                    gridedData.Value[i,j] = gridedData.Value[i,j] + dcoffset
                    if minVal < 0 # also add minvalue to get the lowest value in the dataset to 100
                        gridedData.Value[i,j] = gridedData.Value[i,j] + abs(minVal)
                    end
                end
            end
        end
    end

    # also add the minVal and DC offset to the maxVal so can use it properly later on    
    if minVal < 0
        maxVal += abs(minVal)
    end
    maxVal += dcoffset
    return gridedData, minVal, maxVal, dcoffset
end
