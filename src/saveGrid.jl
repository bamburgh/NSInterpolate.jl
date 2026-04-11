function saveGrid(celldata::CellData, paramd, missingdata; out_file="")
    if out_file == ""
        out_file = paramd["outputFile"]
    end
    ds = NCDataset(out_file, "c")

    value_units = ""
    coord_units = ""
    for mykey in keys(paramd)
        thekey = convert(String, mykey)
        # println("Attribute ", thekey, " = ", paramd[thekey], " (", typeof(paramd[thekey]), ")")
        # println(thekey, " ", paramd[thekey])
        if typeof(paramd[thekey]) in (String , Float64, Int64)
            ds.attrib[thekey] = paramd[thekey]
        end
        if "z_units" == lowercase(mykey)
            value_units = paramd[mykey]
        end
        if "en_units" == lowercase(mykey)
            coord_units = paramd[mykey]
        end
    end

    defDim(ds, "easting", celldata.lenX)
    defDim(ds, "northing", celldata.lenY)

    ncX = defVar(ds, "easting", Float64, ("easting",))
    ncX.attrib["units"] = coord_units
    ncX[:] = celldata.X[:,1]

    ncY = defVar(ds, "northing", Float64, ("northing",))
    ncY.attrib["units"] = coord_units
    ncY[:] = celldata.Y[1,:]

    ncValues = defVar(ds, "Values", Float64, ("easting", "northing"), fillvalue = missingdata)
    ncValues.attrib["units"] = value_units
    ncValues.attrib["field"] = paramd["input_value"]

    myvals = zeros(Float64, celldata.lenX, celldata.lenY)
    for i in 1:celldata.lenX
        for j in 1:celldata.lenY
            if celldata.Flag[i,j] < 0
                myvals[i,j] = missingdata
            else
                myvals[i,j] = celldata.Value[i,j]
            end
        end
    end

    ncValues[:] = myvals

    close(ds)
end
