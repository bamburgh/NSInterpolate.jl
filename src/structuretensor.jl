# now calculate the structure tensor and get the eigenvector for directionality information
function structuretensor(gridedDataDerivIterate, smoothGrid)
    lengthX = gridedDataDerivIterate.lenX
    lengthY = gridedDataDerivIterate.lenY
    m1 = zeros(Float64, lengthX, lengthY)
    v11 = zeros(Float64, lengthX, lengthY)
    v12 = zeros(Float64, lengthX, lengthY)

    for i in 1:lengthX
        for j in 1:lengthY
            # If the data is outside of the grid, then don't use it
            if (gridedDataDerivIterate.Flag[i,j] != -1)
                Ix = 0.0
                if i == 1
                    Ix = smoothGrid[i+1,j] - smoothGrid[i,j]
                elseif (i == lengthX)
                    Ix = smoothGrid[i,j] - smoothGrid[i-1,j]
                else
                    Ix = 0.5 * (smoothGrid[i+1,j] - smoothGrid[i-1,j])
                end
                Iy = 0.0
                if j == 1
                    Iy = smoothGrid[i,j + 1] - smoothGrid[i,j]
                elseif (j == lengthY)
                    Iy = smoothGrid[i,j] - smoothGrid[i,j - 1]
                else
                    Iy = 0.5 * (smoothGrid[i,j + 1] - smoothGrid[i,j - 1])
                end

                # create the structure matrix
                S = zeros(Float64, 2, 2)
                S[1,1] = Ix * Ix
                S[1,2] = Ix * Iy
                S[2,1] = Iy * Ix
                S[2,2] = Iy * Iy

                # calculate eigenvalue for trend determination
                m1[i,j] = 0.5 * (S[1,2] + S[2,2] + sqrt(((S[1,1] - S[2,2]) ^ 2) + (4 * (S[1,2] ^ 2))))

                if (S[1,2] > 0.0)
                    if (S[1,1] == S[2,2])
                        # theta = 1/4 pi
                        theta = 0.25 * pi
                        v12[i,j] = cos(theta)
                        v11[i,j] = -1 * sin(theta)
                    else
                        # 0 < theta < 1/2 pi
                        theta = (atan(2 * S[1,2] / (S[1,1] - S[2,2])) / 2);
                        if (S[1,1] < S[2,2])
                            v11[i,j] = -1 * abs(cos(theta))
                            v12[i,j] = abs(sin(theta))
                        else
                            v12[i,j] = -1 * abs(cos(theta))
                            v11[i,j] = abs(sin(theta))
                        end
                    end
                elseif S[1,2] == 0
                    if S[1,1] > S[2,2]
                        # theta = 0
                        v12[i,j] = 1
                        v11[i,j] = 0
                    else
                        # theta = 1/2 pi
                        v12[i,j] = 0
                        v11[i,j] = 1
                    end
                else
                    if S[1,1] == S[2,2]
                        # theta = 3/4 pi
                        theta = (3 / 4) * pi
                        v12[i,j] = abs(cos(theta))
                        v11[i,j] = abs(sin(theta))
                    else
                        # 1/2 pi < theta < pi
                        theta = atan(2 * S[1,2] / (S[1,1] - S[2,2])) / 2
                        if S[1,1] < S[2,2]
                            v11[i,j] = -1 * abs(cos(theta))
                            v12[i,j] = -1 * abs(sin(theta))
                        else
                            v12[i,j] = -1 * abs(cos(theta))
                            v11[i,j] = -1 * abs(sin(theta))
                        end
                    end
                end
            end
        end
    end
    return m1, v11, v12
end
