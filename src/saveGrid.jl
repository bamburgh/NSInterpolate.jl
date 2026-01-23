function saveGrid(celldata::CellData, paramd, missingdata; out_file="")
    if out_file == ""
        out_file = paramd["outputFile"]
    end
    ds = NCDataset(out_file, "c")

    for mykey in keys(paramd)
        thekey = convert(String, mykey)
        println(thekey, " ", paramd[thekey])
        ds.attrib[thekey] = paramd[thekey]
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
    ncValues.attrib["units"] = "um/s/s"
    ncValues.attrib["field"] = "bouguer gravity (2.67)"

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
