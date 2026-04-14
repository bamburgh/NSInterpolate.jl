"""
 NSInterpolate interpolates data irregularly sampled in 2 dimensions onto a regular grid. It is
 designed originally for aeromagnetic data but is also very effective for airborne gravity and
 airborne gravity gradient data. These data are all sampled along flight-lines that are close to
 being straight lines and it is likely that the method will work well for any data sampled in this
 way.

 The code is a conversion from `c#`, published by T. Naprstek and R. S. Smith (2019, A new method
 for interpolating linear features in aeromagnetic data. Geophysics, 84(3), JM15-JM24).

 I have used version NSI_v8 from github:
     https://github.com/TomasNaprstek/Naprstek-Smith-Interpolation.

 Mark Dransfield 2024
"""
module NSInterpolate

export NSinterp

using Statistics
using Printf
using NCDatasets
using Dates
using JSON

using DimensionalData
using DimensionalData: @dim, XDim, YDim, TimeDim
@dim Fid TimeDim "fiducial"
@dim East XDim "easting"
@dim North YDim "northing"

"""
    NSinterp(param_file::String; verbose=false)

    NSinterp(paramd::Dict; verbose=false)

    NSinterp(observed_data::DimStack, param_file::String; verbose=false)

    NSinterp(observed_data::DimStack, paramd::Dict; verbose=false)

    NSinterp(;input_file::String, east::String, north::String, z::String,
        z_units::String, en_units::String, diff_z::Bool, datum::String, projection::String,
        outputwritebool::Bool, outfile::String, cellSize::Float64, interpDist::Float64,
        maxLoop::Int64, searchStepSize::Float64, cellSizeF::Float64, trendM::Float64,
        autoStop::Bool, angleSearch::Float64, multiSmooth::Float64, spatialSmooth::Bool,
        realGridLocations::Bool, verbose=false
    )

    NSinterp(observed_data::DimensionalData.DimStack; datum::String,
        projection::String, outfile::String, cellSize::Float64, interpDist::Float64,
        maxLoop::Int64, searchStepSize::Float64, cellSizeF::Float64, trendM::Float64,
        autoStop::Bool, angleSearch::Float64, multiSmooth::Float64, spatialSmooth::Bool,
        outputwritebool::Bool, realGridLocations::Bool, verbose=false
    )

 Use control parameters to perform gridding of observed survey data and manage resultant
 grid.
  
 # Arguments
 * `paramd`: A dictionary containing the input parameters (see Notes).
 * `param_file`: The name of the JSON file containing the input parameters.
 * `verbose`: A flag to indicate verbose reporting of progress or not (default).

 # Input Survey Data
 * usually read from the data file named in the input parameters but can be
    input directly as a DimStack of arrays `x`, `y`, and `z` dimensioned by
    time (fiducial).

 # Returns
 * `final_grid`: 

 # Notes
 The parameter file contains values for the following input parameters:

 - input_file: String

    path of data file containing the observed data, either Geosoft XYZ or
    NetCDF4 nc format

 - input_east: String

    the name of the channel in `input_file` containing the `x` values,
    usually eastings.
        
 - input_north: String

    the name of the channel in `input_file` containing the `y` values,
    usually northings.
        
 - input_value: String

    the name of the channel in `input_file` containing the values to grid.
        
 - z_units: String ("")

    the units of the `input_value` data; not used, simply written to output.
        
 - en_units: String ("")

    the units of the `input_east` and `input_north` data for writing to output.

 - diff_z: Bool (false)

    if true, then the `input_value` data are differenced along flight-line before
    gridding.

 - datum: String ("unknown")

    the geographic datum (e.g. WGS84) for the input data; not used, simply
    written to output grid.
        
 - projection: String ("unknown")

    the geographic projection for the input data (e.g. "NUTM17"); not used,
    simply written to output grid.
        
 - outputwritebool: Boolean (true)

    if `false`, no output file is written; if `true`, writes the output
    grid to `outputFile`.
        
 - outputFile: String

    the name of the netCDF4 file to which the grid will be written.
        
 - cellSize: Float

    edge size of each square cell in metres. This is the size of the grid cells
    during the interpolation. N&S generally recommend that the interpolation size be
    0.5 x `cellSizeF`, as this can help smooth out the data. However, this can also
    lead to weak lineaments not trending all the way between two flight lines, as
    they now have "further" to trend. N&S recommend experimentation, but in general
    this value should be half the size, or at most, the same as the output cell size.
        
 - interpDist: Float

    metres away that will be interpolated. Essentially, this determines how
    "far" the method will search for real data when completing the normalisation step.
    In general N&S recommend setting this to 75-100% of the flight line spacing. However,
    if a trend is at a highly acute angle to the flight lines, then it would require a
    much larger interpolation distance.
        
 - maxLoop: Integer (20)

    the number of times the interpolation loop will be processed or maximum number
    of iterations. The maximum number of iterations that you wish the interpolation
    method to go through. This will depend highly on the dataset, but in general
    somewhere between 50 and 100 is enough, particularly when the automatic stopping
    criteria is used.
        
 - searchStepSize: Float

    how much of a cell we will "travel" each search step
        
 - cellSizeF: Float

    resampled final cell size. This is the size of the grid cells at the end
    of the interpolation. This should follow standard interpolation rules, for example
    for aeromagnetic data set to one-quarter to one-fifth the flight line spacing.
        
 - trendM: Float (0.0)

    trend factor. Ranges between 0 and 100, affecting how strongly a lineament
    will be trended (100 is maximum, 0 is minimal trending). In general N&S recommend
    using 100. Irrelevant because a bug is preventing trending from working.
        
 - autoStop: Boolean (true)

    whether or not to auto stop. Set to `true` if you wish to let
    the method determine when there is little change occurring between iterations. If
    set to `false` it will run for `maxLoop` iterations, as determined above.
        
 - angleSearch: Float (10.0)

    the number of degrees it will move each time when searching away from
    the initial eigenvector. When searching for flight line data in the normalization
    process, it is possible that the interpolation distance (described above) will be
    reached before a flight line data cell is found. If this occurs, then the search
    angle (as determined by the trend direction step) will be varied by this "trend
    angle" amount. In general N&S recommend 5-10 degrees. A lower value will be more
    accurate, however may dramatically increase the computation time. A higher value
    will lower the trending accuracy, but be computational quicker.
        
 - multiSmooth: Float (0.0)

    percent amount of smoothing of the multiplier grid before applying the
    normalization process. Ranges between 0(no smoothing) and 100 (maximum smoothing),
    affecting the uniqueness of the normalization values . In general N&S recommend
    setting it to 0, and to only increase if high-frequency "noise" is appearing in
    the interpolation. The subsampling process of the interpolation cell size and
    output cell size will in general be more effective at reducing any high-frequency
    "noise", however in certain cases a high smoothing value (>75) will also help.
        
 - spatialSmooth: Boolean (true)

    whether or not to use spatial smoothing (should usually be `true`). This option
    should almost always be set `true`, as it controls the first step of the iterative
    interpolation process. If turned off, the local derivatives are removed from the
    first step, effectively reducing it to a more simplistic smoothing operator. However,
    it is left as a variable to the user for two specific cases. First, it should always
    be turned off if using latitude/longitude coordinates rather than Northings/Eastings,
    or in any other case where the cell size is set to be less than 1. Second, it can
    be turned off as an additional control to reduce the linear structure within the
    data. In most datasets there will be very little difference between having the
    option on or off, however, if left off, some of the areas with minimal linear
    structure should result in even less linear structure.
        
 - realGridLocations: Boolean (true)

    if `false`, outputs real data in the equi-distance grid cell
    locations. If `true`, then output the real data cells as an average position of
    all real data within the cell. **probably has no effect!!!**

 # Examples
 Here we create a parameter dictionary, then use it to run NSInterp:

 ```julia-repl

 julia> paramd = Dict([
            ("input_file", input_file),
            ("input_east", east),
            ("input_north", north),
            ("input_value", z),
            ("z_units", z_units),
            ("en_units", en_units),
            ("difference_z", diff_z),
            ("datum", datum),
            ("projection", projection),
            ("outputwritebool", outputwritebool),
            ("outputFile", outfile),
            ("cellSize", cellSize),
            ("interpDist", interpDist),
            ("maxLoop", maxLoop),
            ("searchStepSize", searchStepSize),
            ("cellSizeF", cellSizeF),
            ("trendM", trendM),
            ("autoStop", autoStop),
            ("angleSearch", angleSearch),
            ("multiSmooth", multiSmooth),
            ("spatialSmooth", spatialSmooth),
            ("realGridLocations", realGridLocations)
            ])
 julia>  NSinterp(paramd)
 
 NSinterp
  Julia version by Mark Dransfield after Naprstek and Smith
  Version gamma!
  Tue, 24 Jan 2026 10:37:51
  6 threads.
  Julia Version - 1.12.3

 Accessing XYZ data in mydatafile.xyz.

  Found 141 header records
  Found 230 lines
  Found 43 channels

  Starting anisotropic gridding loop, loop counter:  1 2 3 4 5 6 7 8 9 10

 End
 ```

"""
# NSInterp

missingdata = 1.0E10

include("../src/xyzreader.jl")
include("../src/obs_from_geowhizz.jl")
include("../src/localstructs.jl")
include("../src/cells.jl")
include("../src/plotchecks.jl");
include("../src/initial_assign.jl");
include("../src/initial_average.jl");
include("../src/initial_interpolation.jl");
include("../src/offset_positive.jl");
include("../src/alpha_mean.jl");
include("../src/taylor_estimation.jl");
include("../src/smoothgrid.jl");
include("../src/structuretensor.jl");
include("../src/trendcalc.jl");
include("../src/smoothmultiplier.jl");
include("../src/subsample.jl");
include("../src/saveGrid.jl");
include("../src/anisotropic_grid.jl");


function params_from(param_file::String)
    paramd = Dict()
    try
        f = open(param_file, "r")
        paramd = JSON.parse(f)
        close(f)
    catch e
        println("Failed to read JSON: $e")
        println("    using DEFAULTS")
        root = "/Users/markdransfield/Documents/GitHub/AirGravQC/examples/SourceData/"
        root = "/Volumes/MHD Data2024/ActiveSurveys/202403_Xcal_5024_Blackall/Weekly Deliveries/19-04-24/Located/"
        # input_file = root * "Canobie.xyz"
        input_file = root * "2205173_Blackall_AGG_Preliminary.xyz"
        outfile = "Blackall_sm60.nc"
        east = "Easting"
        north = "Northing"
        z = "gD_2P67"
        z_units = "um/s/s"
        en_units = "m"
        diff_z = false

        paramd = Dict([
                ("input_file", input_file),
                ("input_east", east),
                ("input_north", north),
                ("input_value", z),
                ("z_units", z_units),
                ("en_units", en_units),
                ("difference_z", diff_z),
                ("datum", "unknown"),
                ("projection", "unknown"),
                ("outputwritebool", true),
                ("outputFile", outfile),
                ("cellSize", 500),
                ("interpDist", 1200),
                ("maxLoop", 20),
                ("searchStepSize", 0.25),
                ("cellSizeF", 500.0),
                ("trendM", 0.0),
                ("autoStop", true),
                ("angleSearch", 10.0),
                ("multiSmooth", 0.0),
                ("spatialSmooth", true),
                ("realGridLocations", true)
                ])
    end
    return paramd
end


function params_from(input_file::String, east::String, north::String, z::String,
    z_units::String, en_units::String, diff_z::Bool, datum::String, projection::String,
    outputwritebool::Bool, outfile::String, cellSize::Float64, interpDist::Float64,
    maxLoop::Int64, searchStepSize::Float64, cellSizeF::Float64, trendM::Float64,
    autoStop::Bool, angleSearch::Float64, multiSmooth::Float64, spatialSmooth::Bool,
    realGridLocations::Bool
)
    paramd = Dict([
            ("input_file", input_file),
            ("input_east", east),
            ("input_north", north),
            ("input_value", z),
            ("z_units", z_units),
            ("en_units", en_units),
            ("difference_z", diff_z),
            ("datum", datum),
            ("projection", projection),
            ("outputwritebool", outputwritebool),
            ("outputFile", outfile),
            ("cellSize", cellSize),
            ("interpDist", interpDist),
            ("maxLoop", maxLoop),
            ("searchStepSize", searchStepSize),
            ("cellSizeF", cellSizeF),
            ("trendM", trendM),
            ("autoStop", autoStop),
            ("angleSearch", angleSearch),
            ("multiSmooth", multiSmooth),
            ("spatialSmooth", spatialSmooth),
            ("realGridLocations", realGridLocations)
            ])
     return paramd
end


"""
    NSinterp(;input_file::String, east::String, north::String, z::String,
        z_units::String, en_units::String, diff_z::Bool, datum::String, projection::String,
        outputwritebool::Bool, outfile::String, cellSize::Float64, interpDist::Float64,
        maxLoop::Int64, searchStepSize::Float64, cellSizeF::Float64, trendM::Float64,
        autoStop::Bool, angleSearch::Float64, multiSmooth::Float64, spatialSmooth::Bool,
        realGridLocations::Bool, verbose=false
    )

 Use control parameters to perform gridding of observed survey data and manage resultant
 grid.
  
 # Arguments
 * various (see Notes).
 * `verbose`: A flag to indicate verbose reporting of progress or not (default).

 # Notes
 The parameter file contains values for the following input parameters:

 - input_file: String

    path of data file containing the observed data, either Geosoft XYZ or
    NetCDF4 nc format

 - input_east: String

    the name of the channel in `input_file` containing the `x` values,
    usually eastings.
        
 - input_north: String

    the name of the channel in `input_file` containing the `y` values,
    usually northings.
        
 - input_value: String

    the name of the channel in `input_file` containing the values to grid.
        
 - z_units: String ("")

    the units of the `input_value` data; not used, simply written to output.
        
 - en_units: String ("")

    the units of the `input_east` and `input_north` data for writing to output.

 - diff_z: Bool (false)

    if true, then the `input_value` data are differenced along flight-line before
    gridding.

 - datum: String ("unknown")

    the geographic datum (e.g. WGS84) for the input data; not used, simply
    written to output grid.
        
 - projection: String ("unknown")

    the geographic projection for the input data (e.g. "NUTM17"); not used,
    simply written to output grid.
        
 - outputwritebool: Boolean (true)

    if `false`, no output file is written; if `true`, writes the output
    grid to `outputFile`.
        
 - outputFile: String

    the name of the netCDF4 file to which the grid will be written.
        
 - cellSize: Float

    edge size of each square cell in metres. This is the size of the grid cells
    during the interpolation. N&S generally recommend that the interpolation size be
    0.5 x `cellSizeF`, as this can help smooth out the data. However, this can also
    lead to weak lineaments not trending all the way between two flight lines, as
    they now have "further" to trend. N&S recommend experimentation, but in general
    this value should be half the size, or at most, the same as the output cell size.
        
 - interpDist: Float

    metres away that will be interpolated. Essentially, this determines how
    "far" the method will search for real data when completing the normalisation step.
    In general N&S recommend setting this to 75-100% of the flight line spacing. However,
    if a trend is at a highly acute angle to the flight lines, then it would require a
    much larger interpolation distance.
        
 - maxLoop: Integer (20)

    the number of times the interpolation loop will be processed or maximum number
    of iterations. The maximum number of iterations that you wish the interpolation
    method to go through. This will depend highly on the dataset, but in general
    somewhere between 50 and 100 is enough, particularly when the automatic stopping
    criteria is used.
        
 - searchStepSize: Float

    how much of a cell we will "travel" each search step
        
 - cellSizeF: Float

    resampled final cell size. This is the size of the grid cells at the end
    of the interpolation. This should follow standard interpolation rules, for example
    for aeromagnetic data set to one-quarter to one-fifth the flight line spacing.
        
 - trendM: Float (0.0)

    trend factor. Ranges between 0 and 100, affecting how strongly a lineament
    will be trended (100 is maximum, 0 is minimal trending). In general N&S recommend
    using 100. Irrelevant because a bug is preventing trending from working.
        
 - autoStop: Boolean (true)

    whether or not to auto stop. Set to `true` if you wish to let
    the method determine when there is little change occurring between iterations. If
    set to `false` it will run for `maxLoop` iterations, as determined above.
        
 - angleSearch: Float (10.0)

    the number of degrees it will move each time when searching away from
    the initial eigenvector. When searching for flight line data in the normalization
    process, it is possible that the interpolation distance (described above) will be
    reached before a flight line data cell is found. If this occurs, then the search
    angle (as determined by the trend direction step) will be varied by this "trend
    angle" amount. In general N&S recommend 5-10 degrees. A lower value will be more
    accurate, however may dramatically increase the computation time. A higher value
    will lower the trending accuracy, but be computational quicker.
        
 - multiSmooth: Float (0.0)

    percent amount of smoothing of the multiplier grid before applying the
    normalization process. Ranges between 0(no smoothing) and 100 (maximum smoothing),
    affecting the uniqueness of the normalization values . In general N&S recommend
    setting it to 0, and to only increase if high-frequency "noise" is appearing in
    the interpolation. The subsampling process of the interpolation cell size and
    output cell size will in general be more effective at reducing any high-frequency
    "noise", however in certain cases a high smoothing value (>75) will also help.
        
 - spatialSmooth: Boolean (true)

    whether or not to use spatial smoothing (should usually be `true`). This option
    should almost always be set `true`, as it controls the first step of the iterative
    interpolation process. If turned off, the local derivatives are removed from the
    first step, effectively reducing it to a more simplistic smoothing operator. However,
    it is left as a variable to the user for two specific cases. First, it should always
    be turned off if using latitude/longitude coordinates rather than Northings/Eastings,
    or in any other case where the cell size is set to be less than 1. Second, it can
    be turned off as an additional control to reduce the linear structure within the
    data. In most datasets there will be very little difference between having the
    option on or off, however, if left off, some of the areas with minimal linear
    structure should result in even less linear structure.
        
 - realGridLocations: Boolean (true)

    if `false`, outputs real data in the equi-distance grid cell
    locations. If `true`, then output the real data cells as an average position of
    all real data within the cell. **probably has no effect!!!**

"""
function NSinterp(;input_file::String, east::String, north::String, z::String,
    z_units::String, en_units::String, diff_z::Bool, datum::String, projection::String,
    outputwritebool::Bool, outfile::String, cellSize::Float64, interpDist::Float64,
    maxLoop::Int64, searchStepSize::Float64, cellSizeF::Float64, trendM::Float64,
    autoStop::Bool, angleSearch::Float64, multiSmooth::Float64, spatialSmooth::Bool,
    realGridLocations::Bool, verbose=false
)
    paramd = params_from(input_file, east, north, z, z_units, en_units, diff_z, datum,
        projection, outputwritebool, outfile, cellSize, interpDist, maxLoop, searchStepSize,
        cellSizeF, trendM, autoStop, angleSearch, multiSmooth, spatialSmooth, realGridLocations
        )
    return NSinterp(paramd, verbose=verbose)
end


"""
    NSinterp(param_file::String; verbose=false)

 Use control parameters to perform gridding of observed survey data and manage resultant
 grid.
  
 # Arguments
 * `param_file`: The name of the JSON file containing the input parameters.
 * `verbose`: A flag to indicate verbose reporting of progress or not (default).

"""
function NSinterp(param_file::String; verbose=false)
    paramd = params_from(param_file)
    return NSinterp(paramd, verbose=verbose)
end


"""
    NSinterp(paramd::Dict; verbose=false)

 Use control parameters to perform gridding of observed survey data and manage resultant
 grid.
  
 # Arguments
 * `paramd`: A dictionary containing the input parameters (see Notes).
 * `verbose`: A flag to indicate verbose reporting of progress or not (default).

"""
function NSinterp(paramd::Dict; verbose=false)

    println("NSinterp")
    println("  Julia version by Mark Dransfield after Naprstek and Smith")
    println("  Version gamma!")
    println("  ", Dates.format(now(), "e, dd u yyyy HH:MM:SS"))
    println("  ", Threads.nthreads(), " threads.")
    println("  ", "Julia Version - ", VERSION)
    println()

    # Get the observed data
    if occursin(uppercase(split(paramd["input_file"], ".")[end]), "XYZ")
        obs = obs_from_geoxyz(paramd["input_file"];
            n_chan=paramd["input_north"], 
            e_chan=paramd["input_east"], 
            z_chan=paramd["input_value"], 
            outsample=1, 
            verbose=verbose
            )
    elseif  occursin(uppercase(split(paramd["input_file"], ".")[end]), "NC")
        obs = obs_from_geowhizz(paramd["input_file"],
            n_chan=paramd["input_north"], 
            e_chan=paramd["input_east"], 
            z_chan=paramd["input_value"], 
            verbose=verbose)
    else
        suffix = uppercase(split(paramd["input_file"], ".")[end])
        println("error - input data file name must end in either XYZ or NC not $suffix")
        return
    end
    if paramd["difference_z"]
        obs_diff = diff(obs[:down].data[:])
        obs[:down].data[:] = vcat(obs_diff[1], obs_diff)
    end

    return NSinterp(obs, paramd, verbose=verbose)
end


"""
    NSinterp(obs::DimStack, paramd::Dict; verbose=false)

 Use control parameters to perform gridding of observed survey data and manage resultant
 grid.
  
 # Arguments
 * `obs`: A DimStack of geographically located observed data.
 * `paramd`: A dictionary containing the input parameters (see Notes).
 * `verbose`: A flag to indicate verbose reporting of progress or not (default).

"""
function NSinterp(obs::DimensionalData.DimStack, paramd::Dict; verbose=false)

    data = init_xyz(
        true, obs[:east].data, obs[:north].data, obs[:down].data,
        [MapLine(1001., Point(33., 44.), Point(55., 66.))],
        [
        Point(minimum(obs[:east].data), minimum(obs[:north].data)),
        Point(maximum(obs[:east].data), maximum(obs[:north].data))
        ]
        )

    gridedData1, posit, Xmin, Xmax, Ymin, Ymax, X, Y, Value = initial_assign(data, paramd)
    if verbose
        println("\n\nFinished Initial Assignment to Grid")
    end

    println(summary(gridedData1))

    gridedData2 = initial_average(gridedData1, posit)
    if verbose
        println("\n\nFinished Averaging")
        println(summary(gridedData2))
    end

    gridedData3 = initial_interpolation(gridedData2, paramd, Xmin, Ymin, missingdata)
    if verbose
        println("\n\nFinished Initial Interpolating")
        println(summary(gridedData3))
        saveGrid(gridedData3, paramd, missingdata; out_file="init_interpolation.nc")
    end
    
    gridedData4, minVal, maxVal, dcoffset = offset_positive(gridedData3)
    if verbose
        println("\n\nFinished Offsetting to Positive")
    end

    gridedData5 = alpha_mean(gridedData4, paramd)
    if verbose
        println("\n\nFinished Alpha-mean Adjustment")
        println(summary(gridedData5))
        saveGrid(gridedData5, paramd, missingdata; out_file="alpha_mean.nc")
    end

    realReplace = anisotropic_grid(gridedData5, paramd, dcoffset, minVal, maxVal)
    if verbose
        println("\n\nFinished Anisotropic Gridding")
        println(summary(realReplace))
        saveGrid(realReplace, paramd, missingdata; out_file="anisotropic_grid.nc")
    end

    finalData = subsample(X, Y, Xmax, Xmin, Ymax, Ymin, paramd, minVal, Value, realReplace, missingdata, dcoffset)
    if verbose
        println("\n\nFinished Subsampling")
    end

    if paramd["outputwritebool"]
        saveGrid(finalData, paramd, missingdata)
        println("\n\nFinished writing output to $(paramd["outputFile"])")
    end
    println("\n\nNSinterp ended.")

    # return as DimensionalData
    return cell_to_dim(
        finalData,
        name=paramd["input_value"],
        units=paramd["z_units"], 
        datum=paramd["datum"], 
        projection=paramd["projection"]
        )
end

end