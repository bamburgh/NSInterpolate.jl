
struct Point
    x::Float64
    y::Float64
end

struct MapLine
    id::Float64
    start::Point
    last::Point
end

struct XYZ
    dataIsGood::Bool
    X::Array{Float64}
    Y::Array{Float64}
    Value::Array{Float64}
    Line::Array{Union{Nothing, MapLine}}
    bounds::Array{Point} # currently bounding rectangle but could make it a polygon (?)
end

function init_xyz(dataIsGood::Bool=false, X::Array{Float64}=[], Y::Array{Float64}=[], Value::Array{Float64}=[], Line::Array{MapLine}=[], bounds::Array{Point}=[])
    if length(bounds) < 2
        print(@sprintf("ERROR - require at least two points in bounds, only got %d.", length(bounds)))
        dataIsGood = false
    end
    return XYZ(dataIsGood, X, Y, Value, Line, bounds)
end

function summary(xyzdata::XYZ)
    if !xyzdata.dataIsGood
        output = "No valid XYZ data."
    else
        output = @sprintf("\nXYZDATA")
        output *= @sprintf("\nX Field count = %d", length(xyzdata.X))
        output *= @sprintf("\n  min = %.1f", minimum(xyzdata.X))
        output *= @sprintf("\n  max = %.1f", maximum(xyzdata.X))
        output *= @sprintf("\nY Field count = %d", length(xyzdata.Y))
        output *= @sprintf("\n  min = %.1f", minimum(xyzdata.Y))
        output *= @sprintf("\n  max = %.1f", maximum(xyzdata.Y))
        output *= @sprintf("\nValue Field count = %d", length(xyzdata.Value))
        output *= @sprintf("\n  min = %.1f", minimum(xyzdata.Value))
        output *= @sprintf("\n  max = %.1f", maximum(xyzdata.Value))
        output *= @sprintf("\nLine Field count = %d", length(xyzdata.Line))
        output *= @sprintf("\nBounds:")
        output *= @sprintf("\n  X %.1f to %.1f", xyzdata.bounds[1].x, xyzdata.bounds[2].x)
        output *= @sprintf("\n  Y %.1f to %.1f", xyzdata.bounds[1].y, xyzdata.bounds[2].y)
    end
    return output
end

struct Position
    X::Array{Float64}
    Y::Array{Float64}
end

function init_position(nx::Int, ny::Int)
    X = zeros(Float64, nx, ny)
    Y = zeros(Float64, nx, ny)
    return Position(X, Y)
end

function init_position(Xint::Array{Float64}, Yint::Array{Float64})
    
    if size(Xint) != size(Yint)
        print("ERROR - input arrays sizes do not match.")
    end
    if size(Xint)[1] < 2 || size(Xint)[2] < 2
        print("WARNING - input arrays are not 2D.")
    end
    X = Xint
    Y = Yint
    return Position(X, Y)
end

function summary(posData::Position)
    output = @sprintf("\nPOSITION")
    output *= @sprintf("\nX Field count = %d, %d", size(posData.X)[1], size(posData.X)[2])
    output *= @sprintf("\n  min = %.1f", minimum(posData.X))
    output *= @sprintf("\n  max = %.1f", maximum(posData.X))
    output *= @sprintf("\n  spacing = %.1f", abs(posData.X[2] - posData.X[1]))
    output *= @sprintf("\nY Field count = %d, %d", size(posData.Y)[1], size(posData.Y)[2])
    output *= @sprintf("\n  min = %.1f", minimum(posData.Y))
    output *= @sprintf("\n  max = %.1f", maximum(posData.Y))
    output *= @sprintf("\n  spacing = %.1f", abs(posData.Y[2] - posData.Y[1]))
    return output
end

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
    # for  i in 1:mycell.lenX
    #     for j in 1:mycell.lenY
    #         newcell.X[i,j] = mycell.X[i,j]
    #         newcell.Y[i,j] = mycell.Y[i,j]
    #         newcell.Value[i,j] = mycell.Value[i,j]
    #         newcell.Flag[i,j] = mycell.Flag[i,j]
    #     end
    # end
    newcell.X .= mycell.X
    newcell.Y .= mycell.Y
    newcell.Value .= mycell.Value
    newcell.Flag .= mycell.Flag
    newcell.lenX = mycell.lenX
    newcell.lenY = mycell.lenY
    return newcell
end


# USER-INPUT VARIABLES
Base.@kwdef mutable struct Parameters
    datum::String = "WGS84"
    projection::String = "NUTM17"
    outputFile::String = outfile
    cellSize::Float64 = 250 # edge size of each square cell in metres
    interpDist::Float64 = 800 # metres away that will be interpolated
    maxLoop::Int = 20 # the number of times the interpolation loop will be processed 75
    searchStepSize::Float64 = 0.25 # how much of a cell we will "travel" each search step
    cellSizeF::Float64 = 250 # resampled final cell size
    trendM::Float64 = 50 # 100 - median % location (so 0 is no trending, 100 is full trending)
    autoStop::Bool = true # a checkbox of whether or not to auto stop
    angleSearch::Float64 = 10 # the number of degrees it will move each time when searching away from the initial eigenvector
    multiSmooth::Float64 = 60.0 # smooth the multiplier grid before applying the normalization process (0 is no smoothing, 100 is max smoothing) (%)
    spatialSmooth::Bool = true; #a checkbox of whether or not to use spatial smoothing (in almost all cases, should be used)
    outputwritebool::Bool = true; #if false, no output; [outputs in x y value.] if true, outputs in netCDF4 format.
    realGridLocations::Bool = true; #if false, outputs real data in the equi-distance grid cell locations. If true, then output the real data cells as an average position of all real data within the cell.
end

