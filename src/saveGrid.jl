function saveGrid(celldata::CellData, paramd, missingdata; out_file="")
    if out_file == ""
        out_file = paramd["outputFile"]
    end
    ds = NCDataset(out_file, "c")

    for mykey in keys(paramd)
        thekey = convert(String, mykey)
        # println("Attribute ", thekey, " = ", paramd[thekey], " (", typeof(paramd[thekey]), ")")
        # println(thekey, " ", paramd[thekey])
        if typeof(paramd[thekey]) in (String , Float64, Int64)
            ds.attrib[thekey] = paramd[thekey]
        end
    end

    defDim(ds, "easting", celldata.lenX)
    defDim(ds, "northing", celldata.lenY)

    ncX = defVar(ds, "easting", Float64, ("easting",))
    ncX.attrib["units"] = "metres"
    ncX[:] = celldata.X[:,1]

    ncY = defVar(ds, "northing", Float64, ("northing",))
    ncY.attrib["units"] = "metres"
    ncY[:] = celldata.Y[1,:]

    ncValues = defVar(ds, "Values", Float64, ("easting", "northing"), fillvalue = missingdata)
    if "units" in lowercase(keys(paramd))
        ncValues.attrib["units"] = paramd["units"]
    else
        ncValues.attrib["units"] = ""
    end
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
