"""
	obs_from_geoxyz(xyz_file)

    Read a set of XYZ observation locations from a Geosoft XYZ file.
"""
function obs_from_geoxyz(xyz_file; n_chan="y", e_chan="x", z_chan="z", missing_value=-1.0E-64, outsample=1, verbose=false)
    # fspec = FormatSpec("E")
    missing_str = "-1.0E-64"#fmt(fspec, missing_value)
    # count number of header records, number of flight lines in xyz_file and number of channels
    (num_head_recs, num_lines, num_channels, field_precisions) = xyz_count(xyz_file, verbose)
    
    # get channel names ...
    channelnames = xyz_channels(xyz_file, num_head_recs, num_channels; verbose=verbose)

    # ... and check we have the ones we want.
    # We expect northing, easting and vertical positions. We need n, e, d (d = - vertical).
    n_index = findfirst(isequal(n_chan), channelnames)
    if isnothing(n_index)
        println(@sprintf("Northing channel %s not found. Channel names were:", n_chan))
        println(channelnames)
        return [0. 0. 0.]
    end
    e_index = findfirst(isequal(e_chan), channelnames)
    if isnothing(e_index)
        println(@sprintf("Easting channel %s not found. Channel names were:", e_chan))
        println(channelnames)
        return [0. 0. 0.]
    end
    z_index = findfirst(isequal(z_chan), channelnames)
    if isnothing(z_index)
        println(@sprintf("Height channel %s not found. Channel names were:", z_chan))
        println(channelnames)
        return [0. 0. 0.]
    end

    # get line numbers and lengths
    (line_ids, num_fids) = xyz_lines(xyz_file, num_lines)

    total_fids = sum(num_fids)

    # We will store our observations in a DimStack of 3 Vector DimArrays with names and units.

    my_fids = collect(1:1:total_fids)
    x1 = DimArray(zeros(total_fids), (Fid(my_fids),), name="north", metadata=Dict("units"=>"m", "orig_name"=>n_chan))
    x2 = DimArray(zeros(total_fids), (Fid(my_fids),), name="east", metadata=Dict("units"=>"m", "orig_name"=>e_chan))
    x3 = DimArray(zeros(total_fids), (Fid(my_fids),), name="down", metadata=Dict("units"=>"m", "orig_name"=>z_chan))
    obs = DimStack(x1, x2, x3)

    # obs = Obs(
    #     n_chan, e_chan, z_chan,
    #     "m", "m", "m",
    #     zeros(total_fids), zeros(total_fids), zeros(total_fids)
    #     )

    open(xyz_file) do xyzf
        fids_read = 0
        datarecs_read = 0
        for xyz_rec in eachline(xyzf)
            if startswith(xyz_rec, "/")
            elseif start_upper(xyz_rec, "LINE") || start_upper(xyz_rec, "TIE")
            else
                # check & clean xyz_rec
                rec_str = split(xyz_rec)
                if length(rec_str) != num_channels
                    println("ERROR - Line $line_number - Number of channel names mis-match to size of data ")
                    println("Record has $(length(rec_str)) elements. Number of channels = $num_channels")
                    println("Error occurred after $recs_of_line_read records in line.")
                    break
                end
                datarecs_read += 1
                if rem(datarecs_read, outsample) == 0
                    for ii = 1:num_channels
                        if rec_str[ii] == "*"
                            rec_str[ii] = missing_str
                        end
                    end
                    rec_values = parse.([Float64], rec_str)
                    fids_read += 1
                    obs[:north].data[fids_read] = rec_values[n_index]
                    obs[:east].data[fids_read] = rec_values[e_index]
                    obs[:down].data[fids_read] = rec_values[z_index]
                end
            end
        end
    end

    return obs
end

"""
    xyz_count(file_name::String)

    Count the number of header records, flight lines and channels in XYZ file.

    Also the decimal precision of each channel in the Geosoft XYZ file.

    params
    ----------
    file_name ::String
        The name of the Geosoft XYZ file.

    Returns
    -------
    num_head_recs ::Int
        The .
    num_lines ::Int
        The .
    num_channels ::Int
        The .
    field_precisions ::Int
        The .
"""
function xyz_count(file_name::String, verbose=false)
    num_lines = 0
    num_head_recs = 0
    num_channels = 0
    need_first_data_rec = true
    field_precisions = Int64[]
    xyz_filename = splitpath(file_name)[end]

    open(file_name) do fid
        println(@sprintf("Accessing XYZ data in %s.", xyz_filename))
        if verbose
            println("\n\nFirst few records are:")
        end
        for file_line in eachline(fid)
            if num_lines == 0 && verbose
                # Always useful to see the first few records in the file.
                println(file_line)
            end
            if file_line[1] == '/'
                num_head_recs += 1
            elseif start_upper(file_line, "LINE") || start_upper(file_line, "TIE")
                num_lines += 1
            # Need every value in the row to be a decimal number (no '*' dummies)
            elseif occursin("*", file_line)
            elseif num_lines > 0 && need_first_data_rec
                first_data_rec = split(file_line)
                num_channels = length(first_data_rec)
                for ii = 1:length(first_data_rec)
                    bits = split(first_data_rec[ii], '.')
                    if bits[1] == first_data_rec[ii]
                        push!(field_precisions, 0)
                    else
                        push!(field_precisions,length(bits[2]))
                    end
                end
                need_first_data_rec = false
            else
            end
        end
    end
    println("\n  Found $num_head_recs header records")
    println("  Found $num_lines lines")
    println("  Found $num_channels channels\n")
    return (num_head_recs, num_lines, num_channels, field_precisions)
end

"""
    xyz_lines(file_name::String, num_lines)

    Returns the line numbers (`line_ids`), and the number of fiducials in each line.
    A helper funtion for `xyz_to_whizz`.

    params
    ----------
    whizz_file ::String
        The name of the geoWhizz file.

    Returns
    -------
    None
"""
function xyz_lines(file_name::String, num_lines)

    # initialise counter for the number of fids per line and storage for lines, flights and dates
    num_fids = zeros(Int64, num_lines)
    line_ids = String[]
    line_ctr = 1
    # flight_nos = zeros(Int64, num_lines)
    # flight_dates = [[1, 1, 1980] for k in range(num_lines)]
    # flight_no = 0
    # flight_date = [1, 1, 1980]

    open(file_name) do fid
        for file_line in eachline(fid)
            if length(file_line) < 4
            elseif file_line[1] == '/'
                if file_line[2] != '/'  # already got channel names
                elseif start_upper(file_line, "//FLIGHT")
                    flight_no = parse(Int64, split(file_line)[2])
                elseif start_upper(file_line, "//DATE")
                    test = split(file_line)[2]
                    y, m, d = split(test, '/')
                    flight_date = [parse(Int64, d), parse(Int64, m), parse(Int64, y)]
                end
            elseif start_upper(file_line, "LINE") || start_upper(file_line, "TIE")
                push!(line_ids, split(file_line)[2])
                # flight_nos[line_ctr] = flight_no
                # flight_dates[line_ctr] = flight_date
                line_ctr += 1
            # elseif file_line.count('*') != 0
            #     continue  # print('Skip record for dummy')
            else
                num_fids[line_ctr-1] += 1
            end
        end
    end
    return line_ids, num_fids
end

"""
    xyz_channels(file_name::String, num_head_recs::Integer, num_channels::Integer)

 Get the names of the channels in a Geosoft XYZ `file_name`. The algorithm checks
 `num_head_recs` header records. If it finds one with a number of words equal to
 the number of channels, `num_channels`, then it returns those words as an array
 of channel names.
"""
function xyz_channels(file_name::String, num_head_recs::Integer, num_channels::Integer; verbose=false)

    channelnames = repeat([""],num_channels)
    header_rec = 0
    open(file_name) do fileid
        for file_line in eachline(fileid)
            temp_names = split(lstrip(file_line, '/'))
            header_rec += 1
            if header_rec > num_head_recs
                channelnames = [""]
                println("Error - can't find header record with $num_channels channel names.")
                break
            elseif length(temp_names) == num_channels
                channelnames = temp_names
                if verbose
                    for ii = 1:length(channelnames)
                        println("$(channelnames[ii])")
                    end
                end
                break
            end
        end
    end
    return channelnames
end

start_upper(my_string::String, start_string::String) = startswith(uppercase(lstrip(my_string)), start_string)


