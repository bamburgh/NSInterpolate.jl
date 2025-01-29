function findbase(mymatrix, missing)
    mx = copy(mymatrix)
    if mean(mx) > missing
        mx[mymatrix .== missing] .= 1000000000
        base = minimum(mx)
    else
        mx[mymatrix .== missing] .= -1000000000
        base = maximum(mx)
    end
    return base
end

function plotmatrix(mymatrix, missing, mytitle)
    mx = copy(mymatrix)
    mx[mymatrix .== missing] .= findbase(mymatrix, missing)
    return heatmap(transpose(mx), title=mytitle, tickfontsize=6, titlefontsize=9, transpose=true, framestyle=:box, aspect_ratio=:equal)
end

function cleanXYZ(myxyz, missing)
    mx = myxyz
    mx.Value[myxyz.Value .== missing] .= findbase(myxyz.Value, missing)
    return mx
end
    

function plotXYZ(myxyz, missing, truecoords)
    px = plotmatrix(myxyz.X, missing, "X")
    py = plotmatrix(myxyz.Y, missing, "Y")
    if truecoords
        mx = cleanXYZ(myxyz, missing)
        pv = heatmap(mx.X[:,1], mx.Y[1,:], transpose(mx.Value), title="Value", tickfontsize=6, titlefontsize=9, framestyle=:box, aspect_ratio=:equal)
    else
        pv = plotmatrix(myxyz.Value, missing, "Value")
    end
    pf = plotmatrix(myxyz.Flag, -2, "Flag")

    plot(px, py, pf, pv, layout=(2,2), size=(600, 800))
end
