"""
The code is a conversion from `c#` of that published by T. Naprstek and R. S. Smith (2019, 
A new method for interpolating linear features in aeromagnetic data. Geophysics, 84(3), 
JM15-JM24). I have used version NSI_v8 from github:
    https://github.com/TomasNaprstek/Naprstek-Smith-Interpolation.
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

missingdata = 1.0E10

include("../src/xyzreader.jl")
include("../src/localstructs.jl")
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

function read_json(file::String)::Dict
    open(file, "r") do f
        return JSON.parse(f)
    end
end

function NSinterp(;xyz_file, east, north, z, datum, projection, outfile, cellSize, interpDist,
	maxLoop, searchStepSize, cellSizeF, trendM, autoStop, angleSearch, multiSmooth, spatialSmooth,
	outputwritebool, realGridLocations)
	paramd = Dict([
	        ("input_xyz_file", xyz_file),
	        ("input_xyz_east", east),
	        ("input_xyz_north", north),
	        ("input_xyz_value", z),
	        ("datum", datum),
	        ("projection", projection),
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
	        ("outputwritebool", outputwritebool),
	        ("realGridLocations", realGridLocations)
	        ])
	NSinterp(paramd)
end

"""
    NSinterp(param_file)

Read parameters from `param_file` into a dictionary, `paramd`, and call `NSinterp(paramd)`

A more detailed explanation can go here, perhaps a ref to the paper.

# Arguments
* `param_file`: The name of the JSON file containing the input parameters.

# Notes
The parameter file contains values for the following input parameters:

    input_xyz_file: Geosoft XYZ file containing the observed data
    input_xyz_east: the name of the channel in `input_xyz_file` containing the eastings
    input_xyz_north: the name of the channel in `input_xyz_file` containing the northings
	input_xyz_value: the name of the channel in `input_xyz_file` containing the values to grid
    datum: the geographic datum (e.g. WGS84) for the input data
    projection: the geographic projection for the input data (e.g. "NUTM17")
    outputFile: the name of the netCDF4 file that the grid will be written to
	cellSize: edge size of each square cell in metres
	interpDist: metres away that will be interpolated
	maxLoop: the number of times the interpolation loop will be processed
	searchStepSize: how much of a cell we will "travel" each search step
	cellSizeF: resampled final cell size
	trendM: 100 - median % location (so 0 is no trending, 100 is full trending)
	autoStop: a checkbox of whether or not to auto stop
	angleSearch: the number of degrees it will move each time when searching away from the initial eigenvector
	multiSmooth: smooth the multiplier grid before applying the normalization process (0 is no smoothing, 100 is max smoothing) (%)
	spatialSmooth: a checkbox of whether or not to use spatial smoothing (in almost all cases, should be used)
	outputwritebool: if 0, outputs in x y value. if 1, outputs in a format easy for importing into Oasis Montaj.
	realGridLocations: if 0, outputs real data in the equi-distance grid cell locations. If 1, then output the real data cells as an average position of all real data within the cell.

# Examples
```julia
julia>  NSinterp("tokens.json")

NSinterp
  Julia version by Mark Dransfield after Naprstek and Smith
  Version gamma!
  Tue, 28 Jan 2025 10:37:51
  6 threads.
  Julia Version - 1.8.2

Accessing XYZ data in 2205173_Blackall_AGG_Preliminary.xyz.

  Found 141 header records
  Found 230 lines
  Found 43 channels

  Starting anisotropic gridding loop, loop counter:  1 2 3 4 5 6 7 8 9 10

End
```
"""
function NSinterp(param_file::String)
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
		# xyz_file = root * "Canobie.xyz"
		xyz_file = root * "2205173_Blackall_AGG_Preliminary.xyz"
		outfile = "Blackall_sm60.nc"
		east = "Easting"
		north = "Northing"
		z = "gD_2P67"

		paramd = Dict([
		        ("input_xyz_file", xyz_file),
		        ("input_xyz_east", east),
		        ("input_xyz_north", north),
		        ("input_xyz_value", z),
		        ("datum", "unknown"),
		        ("projection", "unknown"),
		        ("outputFile", outfile),
		        ("cellSize", 500),
		        ("interpDist", 1200),
		        ("maxLoop", 20),
		        ("searchStepSize", 0.25),
		        ("cellSizeF", 500),
		        ("trendM", 50),
		        ("autoStop", true),
		        ("angleSearch", 10),
		        ("multiSmooth", 100.0),
		        ("spatialSmooth", true),
		        ("outputwritebool", true),
		        ("realGridLocations", true)
		        ])
	end
	NSinterp(paramd)
end

function NSinterp(paramd::Dict; verbose=false)

	println("NSinterp")
	println("  Julia version by Mark Dransfield after Naprstek and Smith")
	println("  Version gamma!")
	println("  ", Dates.format(now(), "e, dd u yyyy HH:MM:SS"))
	println("  ", Threads.nthreads(), " threads.")
	println("  ", "Julia Version - ", VERSION)
	println()

	obs = obs_from_geoxyz(paramd["input_xyz_file"]; 
		n_chan=paramd["input_xyz_north"], 
		e_chan=paramd["input_xyz_east"], 
		z_chan=paramd["input_xyz_value"], 
		outsample=1, 
		verbose=false
		)

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

	gridedData2 = initial_average(gridedData1, posit)
    if verbose
    	println("\n\nFinished Averaging")
    end

	gridedData3 = initial_interpolation(gridedData2, paramd, Xmin, Ymin, missingdata)
    if verbose
    	println("\n\nFinished Initial Interpolating")
    end

	gridedData4, minVal, maxVal, dcoffset = offset_positive(gridedData3)
    if verbose
    	println("\n\nFinished Offsetting to Positive")
    end

	gridedData5 = alpha_mean(gridedData4, paramd)
    if verbose
    	println("\n\nFinished Alpha-mean Adjustment")
    end

	realReplace = anisotropic_grid(gridedData5, paramd, dcoffset, minVal, maxVal)
    if verbose
    	println("\n\nFinished Anisotropic Gridding")
    end

	finalData = subsample(X, Y, Xmax, Xmin, Ymax, Ymin, paramd, minVal, Value, realReplace, missingdata, dcoffset)
    if verbose
    	println("\n\nFinished Subsampling")
    end

	if paramd["outputwritebool"]
	    saveGrid(finalData, paramd, missingdata)
	    if verbose
	    	println("\n\nFinished writing output")
	    end
	end
	println("NSinterp ended.")
end

end