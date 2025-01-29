"""
The code is a conversion from `c#` of that published by T. Naprstek and R. S. Smith (2019, 
A new method for interpolating linear features in aeromagnetic data. Geophysics, 84(3), 
JM15-JM24). I have used version NSI_v8 from github:
    https://github.com/TomasNaprstek/Naprstek-Smith-Interpolation.
"""
Module NSinterpolate

export NSinterp
# using Plots
using Statistics
# define fiducial dimension
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

function NSinterp(param_file::String)
	# read parameters from `param_file` into paramd and call NSinterp(paramd::Dict)
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