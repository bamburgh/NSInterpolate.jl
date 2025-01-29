"""
Check through all grid cells, and if a cell has more than one reading in it:
average the value over the number of readings. Also average the relative positions
and restore the offset, so that the mean value is at the mean position and not at
the centre, or a corner of, the cell. The flag for the cell is then set to 1.
"""
function initial_average(gridedData, posit)
    for i in 1:gridedData.lenX
        for j in 1:gridedData.lenY
            if (gridedData.Flag[i, j] >= 1)
                gridedData.X[i, j] = gridedData.X[i, j] / (gridedData.Flag[i, j]) + posit.X[i, j]
                gridedData.Y[i, j] = gridedData.Y[i, j] / (gridedData.Flag[i, j]) + posit.Y[i, j]
                gridedData.Value[i, j] = gridedData.Value[i,j] / (gridedData.Flag[i, j])
                # Set it back to 1, showing it is populated by a "true" value
                gridedData.Flag[i, j] = 1
            end
        end
    end
    return gridedData
end
