

mutable struct CellData
    # All 2-D arrays, e.g. [3][7]
    #  Flag = -1 => cell is too far from data to be used
    #       =  0 => no data n cell, expect interpolator to fill
    #       =  1 => has one 'true' value, either from one measured, or by averaging several in cell
    #       >  0 => number is number of data points in cell
    X::Matrix{Float64}
    Y::Matrix{Float64}
    Value::Matrix{Float64}
    Flag::Matrix{Int}
        
    # 
    #  0 = Left
    #  1 = Up left
    #  2 = Up
    #  3 = Up right
    #  4 = Right
    #  5 = Down right
    #  6 = Down
    #  7 = Down left
    # 
    # distance to real data cell in 8 directions, close[x,y,0]=Left,distance
    close::Array{Int}
    # x,y location of 1 cell away data in each direction (i.e. already taken into account edges, etc.), close1[x,y,3,0]=UR,x, close1[x,y,3,1]=UR,y
    close1::Array{Int}
    
    lenX::Int
    lenY::Int
end


function init_cell(nx::Int, ny::Int)
    X = zeros(Float64, nx, ny)
    Y = zeros(Float64, nx, ny)
    Value = zeros(Float64, nx, ny)
    Flag = zeros(Int, nx, ny)
    close = zeros(Int, nx, ny, 8)
    close1 = zeros(Int, nx, ny, 8, 2)
    return CellData(X, Y, Value, Flag, close, close1, nx, ny)
end


function init_cell(Xint::Array{Float64}, Yint::Array{Float64}, Valueint::Array{Float64}, Flagint::Array{Int}, close::Array{Int}, close1::Array{Int})
    X = Xint
    Y = Yint
    Value = Valueint
    Flag = Flagint
    close = close
    close1 = close1
    if abs(X[2,1] - X[1,1]) > abs(X[1,2] - X[1,1])
        nx = size(celldata.X)[1]
        ny = size(celldata.Y)[2]
    else
        nx = size(celldata.X)[2]
        ny = size(celldata.Y)[1]
    end
    return CellData(X, Y, Value, Flag, close, close1, nx, ny)
end


function summary(celldata::CellData)
    output = @sprintf("\nCELLDATA")
    output *= @sprintf("\nX Field count = %d, %d", size(celldata.X)[1], size(celldata.X)[2])
    output *= @sprintf("\n  min = %.1f", minimum(celldata.X))
    output *= @sprintf("\n  max = %.1f", maximum(celldata.X))
    output *= @sprintf("\nY Field count = %d, %d", size(celldata.Y)[1], size(celldata.Y)[2])
    output *= @sprintf("\n  min = %.1f", minimum(celldata.Y))
    output *= @sprintf("\n  max = %.1f", maximum(celldata.Y))
    output *= @sprintf("\nValue Field count = %d, %d", size(celldata.Value)[1], size(celldata.Value)[2])
    output *= @sprintf("\n  min = %.1f", minimum(celldata.Value))
    output *= @sprintf("\n  max = %.1f", maximum(celldata.Value))
    return output
end


function copy_cell(mycell, newcell)
    newcell.X .= mycell.X
    newcell.Y .= mycell.Y
    newcell.Value .= mycell.Value
    newcell.Flag .= mycell.Flag
    newcell.lenX = mycell.lenX
    newcell.lenY = mycell.lenY
    return newcell
end


"""

"""
function cell_to_dim(mycell; name="", units="", projectname="", datum="", projection="", description="")
    eastings = mycell.X[:,1]
    northings = mycell.Y[1,:]
    newdata = [x == missingdata ? missing : x for x in mycell.Value]
    mydata = DimArray(
        newdata, #mycell.Value,
        (X(eastings), Y(northings));
        name=name#:gravity,
        metadata=Dict(
            "x" => "easting",
            "y" => "northing",
            "Units" => units,
            "Project" => projectname,
            "Datum" => datum,
            "Projection" => projection,
            "Description" => description
    ))
    return mydata
end
