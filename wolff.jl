# Step (RECURSIVE) of Wolff's Algorith
function wolffclusterstep!(spinmatrix::Array{Int, 2},
                    cluster::BitArray{2},
                    i::Int, 
                    j::Int; 
                    temp::Float64 = 1.0, 
                    h::Float64    = 0.0)

    if temp == 0 error("The Temperature can't be ZERO!!") end

    cluster[i,j] = true
    myspin = spinmatrix[i,j]

    for (ni, nj) in neighbors(spinmatrix, i, j)
        if cluster[ni,nj] == false && spinmatrix[ni,nj] == myspin #dunno if the order matters
           if rand() < 1 - exp(-2/temp) wolffclusterstep!(spinmatrix, cluster, ni, nj, temp=temp, h=h) end
        end
    end
end

# Wolff's Algorithm
function wolff!(spinmatrix::Array{Int, 2};
                temp::Float64 = 1.0,
                h::Float64    = 0.0,
                maxit::Int    = 1000,
                plot::Bool    = true,
                verbose::Bool = true)

    # Generating the postion randomically
    x , y = rand(1:size(spinmatrix, 1)) , rand(1:size(spinmatrix, 12))
    cluster = falses(size(spinmatrix))
    mi = Array(Float64, 0)
    xi = 1:maxit

    for iter in xi
        wolffclusterstep!(spinmatrix, cluster, x, y, temp=temp, h=h)
        flip!(spinmatrix, cluster)
        push!(mi, abs(magnetization(spinmatrix)))
        println("$(mi[end])")
    end

    if verbose println("Finished with magnetization $(mi[end])") end
    if plot
        PyPlot.plot(xi, mi, "o", color="blue")
        PyPlot.title("Wolff on Ising for T=$temp")
        PyPlot.xlabel("Number of Iterations")
        PyPlot.ylabel("Magnetization")        
        PyPlot.ylim(0,1.1)
        PyPlot.savefig("Plots/Wolff/wolff_mag_$temp.png")
        PyPlot.close()
    end

    return mi[end]
end
