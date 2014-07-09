# Step (RECURSIVE) of Wolff's Algorith
function wolffstep!(spinmatrix::Array{Int, 2},
                    i::Int, 
                    j::Int; 
                    temp::Float64 = 1.0, 
                    h::Float64    = 0.0,
                    maxit::Int   = 10000)

    if temp == 0 error("The Temperature can't be ZERO!!") end

    # Calculating neighbors' spin sum
    mclose = sum(spinneighbors(spinmatrix, i, j))

    ΔE = 2(h+mclose)*spinmatrix[i,j] # energy necessary to flip the spin
    if ΔE < 0
        flip!(spinmatrix, i, j) # flip spin certainly
    elseif rand() < exp(-ΔE/temp)
        flip!(spinmatrix, i, j) # flip spin with probability
    end

    # Recursive flip neighbor with same spin
    for n in neighbors(spinmatrix, i, j)
        if spinmatrix[n[1], n[2]] == spinmatrix[i,j]
            wolffstep!(spinmatrix, n[1], n[2], temp=temp, h=h, maxit=maxit)
        end
    end
end

# Wolff's Algorithm
function wolff!(spinmatrix::Array{Int, 2};
                temp::Float64 = 1.0,
                h::Float64    = 0.0,
                maxit::Int    = 5000,
                plot::Bool    = true,
                verbose::Bool = true)

    n = size(spinmatrix, 1)
    x , y = rand(1:n) , rand(1:n) # Generating the postion randomically
    mi = Array(Float64, 0)
    xi = 1:maxit

    for iter in xi
        wolffstep!(spinmatrix, x, y, temp=temp, h=h, maxit=maxit)
        push!(mi, abs(magnetization(spinmatrix)))
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
