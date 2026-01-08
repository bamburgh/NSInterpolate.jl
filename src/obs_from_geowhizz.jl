"""
	obs_from_geowhizz(whizz_file)

    Read a set of observations and locations from a geoWhizz file.
"""
function obs_from_geowhizz(whizz_file; n_chan="y", e_chan="x", z_chan="z", verbose=false)

    northing, easting, z_data = get_located_data(whizz_file, z_chan;
        north_chan=n_chan, east_chan=e_chan, lines=[])

    total_fids = length(northing)

    # We will store our observations in a DimStack of 3 Vector DimArrays with names and units.

    my_fids = collect(1:1:total_fids)
    x1 = DimArray(northing, (Fid(my_fids),), name="north", metadata=Dict("units"=>"m", "orig_name"=>n_chan))
    x2 = DimArray(easting, (Fid(my_fids),), name="east", metadata=Dict("units"=>"m", "orig_name"=>e_chan))
    # need to fix the name and units here!!! ToDo
    x3 = DimArray(z_data, (Fid(my_fids),), name="down", metadata=Dict("units"=>"m", "orig_name"=>z_chan))
    obs = DimStack(x1, x2, x3)

    return obs
end


"""
    get_located_data(whizz_file::String, line, z_chan; north_chan="", east_chan"")

 Returns located data for a given survey line.
"""
function get_located_data(whizz_file, z_chan; north_chan="", east_chan="", lines=[])
    northing = Float64[]
    easting = Float64[]
    z_data = Float64[]

    NCDataset(whizz_file) do df
        lines_group = df.group["Lines"]

        if north_chan == ""
            n_name = df.attrib["northing"]
        else
            n_name = north_chan
        end
        if east_chan == ""
            e_name = df.attrib["easting"]
        else
            e_name = east_chan
        end
        if lines == []
            lines = keys(lines_group.group)
        end

        for line in lines
            line_group = lines_group.group[line]
            northing = append!(northing, line_group[n_name][:])
            easting = append!(easting, line_group[e_name][:])
            z_data = append!(z_data, line_group[z_chan][:])
        end
    end
    return northing, easting, z_data
end
