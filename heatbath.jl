# Step of the HeatBath Algorithm
function heatbathstep!(spinmatrix::Array{Int64,2}; temp::Float64 = 1.0, h::Float64=0.0, verbose::Bool=true)
    if temp == 0 error("The Temperature can't be ZERO!!") end

    n = size(spinmatrix, 1)
    x , y = rand(1:n) , rand(1:n) # Generating the position randomically

    # Calculating neighbors' spin sum
    magnear = sum(spinneighbors(spinmatrix, x, y))

    # Calculating the energy on this environment
    energy_up = -magnear - h
    energy_down = magnear + h

    # Calculating the spin new "probabilities"
    prob_up = exp(-energy_up/temp)
    normalization = prob_up + exp(-energy_down/temp)

    # Changing spin with the new probabilities
    spinmatrix[x,y] = rand()*normalization < prob_up ? 1 : -1
    if verbose 
        println("Changed spin ($(int(x)),$(int(y))) to $(spinmatrix[i,j])") 
    end
end

# HeatBath Algorithm for a Given Graph at a Given Temperature
function heatbath!(spinmatrix::Array{Int64,2}; temp::Float64=1.0, h::Float64=0.0, maxit::Int=20000, plot::Bool=true, verbose::Bool=true)
    xi = 1:maxit
    mi = Array(Float64, 0)
    
    for iter in xi
        heatbathstep!(spinmatrix, temp=temp, h=h, verbose=false)
        push!(mi, abs(magnetization(spinmatrix)))
    end
    
    if verbose println("Finished with magnetization $(mi[end])") end
    if plot
        PyPlot.plot(xi, mi, "o", color="blue")
        PyPlot.title("HeatBath on Ising for T=$temp")
        PyPlot.xlabel("Number of Iterations")
        PyPlot.ylabel("Magnetization")   
        PyPlot.savefig("Plots/HeatBath/hb_mag_$temp.png")
        PyPlot.close()
    end

    return mi[end]
end

# HeatBath Algorithm for many Graphs at a Given Temperature
function heattemp(n::Int, temp::Float64; qtd::Int=200, h::Float64=0.0, maxit::Int=20000)
    mi = Array(Float64,0)

    for i in 1:qtd
        ensemble = spingrid(n)
        push!(mi, heatbath!(ensemble, temp=temp, h=h, maxit=maxit, plot=false, verbose=false))
    end
    println("Finished temperature $temp")

    return mean(mi) 
end

# HeatBath Algorithm for many Graphs at several Temperatures
function heatrange(n::Int; maxtemp::Float64=6.0, qtd::Int=200, h::Float64=0.0, maxit::Int=20000)
    mα = Array(Float64,0)

    for i in 0.1:0.1:maxtemp
       push!(mα,heattemp(n, i, qtd=qtd, h=h, maxit=maxit))
    end

    PyPlot.plot(0.1:0.1:maxtemp, mα, "-", color="red")
    PyPlot.title("Magnetization over Temperatures with HeatBath")
    PyPlot.savefig("Plots/HeatBath/heatbath_$(n)grid_$(int(maxtemp))")
    PyPlot.close()

    return mα
end