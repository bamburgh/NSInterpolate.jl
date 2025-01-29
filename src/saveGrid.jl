function saveGrid(celldata::CellData, paramd, missingdata)
    ds = NCDataset(paramd["outputFile"], "c")

    ds.attrib["geographic_datum"] = "WGS84"
    ds.attrib["utm_zone"] = "UTM54"

    for mykey in collect(keys(paramd))
        thekey = convert(String, mykey)
        ds.attrib[thekey] = "test"#paramd[mykey]
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
