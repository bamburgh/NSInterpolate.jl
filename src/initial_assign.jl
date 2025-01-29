"""
Constructs a zero grid (CellData), then for each observation, adds its value 
and location, relative to the bottom-left corner of the cell, to the grid,
together with a flag. The flag's value is set as the number of observations
inside the cell.

(Later, one can average the contained observations and positions.)

Returns the grid, the positions of the bottom left corners of the grid,
and the positions (with their extrema) and values of the observations.
"""

function initial_assign(data::XYZ, params::Dict)
    # working copy
    X = data.X
    Y = data.Y
    Value = data.Value
    datalength = length(X)

    # create grid arrays
    Xmax = maximum(X)
    Xmin = minimum(X)
    Ymax = maximum(Y)
    Ymin = minimum(Y)

    # ... then calc number of cells in starting grid
    lengthX = Int(ceil(((Xmax - Xmin) / params["cellSize"])))
    lengthY = Int(ceil(((Ymax - Ymin) / params["cellSize"])))
    # # (if the division is perfect then we need to add one)
    if mod((Xmax - Xmin), params["cellSize"]) == 0.0
        lengthX += 1
    end
    if mod((Ymax - Ymin), params["cellSize"]) == 0.0
        lengthY += 1
    end

    # #  ... now create starting arrays to hold starting values
    gridedData = init_cell(lengthX, lengthY)
    # print(summary(gridedData))

    posit = init_position(lengthX, lengthY)
    # print(summary(posit))

    # assign initial values from input data
    for k in 1:datalength
        # Round down Ints so the xPos will be in the grid referenced to bottom left point.
        xPos = Int(floor((X[k] - Xmin) / params["cellSize"])) + 1 # + 1
        yPos = Int(floor((Y[k] - Ymin) / params["cellSize"])) + 1 # + 1
        # find bottom left position
        posit.X[xPos, yPos] = (xPos * params["cellSize"]) + Xmin
        posit.Y[xPos, yPos] = (yPos * params["cellSize"]) + Ymin
        # Add the x & y pos of the reading (Its distance from cell's bottom left.) to the cell. 
        gridedData.X[xPos, yPos] += X[k] - posit.X[xPos, yPos]
        gridedData.Y[xPos, yPos] += Y[k] - posit.Y[xPos, yPos]
        # Add the value of that reading to the cell.
        gridedData.Value[xPos, yPos] += Value[k]
        # Track number of readings assigned to the cell.
        gridedData.Flag[xPos, yPos] += 1
    end
    return gridedData, posit, Xmin, Xmax, Ymin, Ymax, X, Y, Value
end
