# Step of the HeatBath Algorithm
function heatbathstep!(spinmatrix::Array{Int64,2}; temp::Float64 = 1.0, h::Float64=0.0)
    if temp == 0 error("The Temperature can't be ZERO!!") end

    # Inicializations
    n = size(spinmatrix, 1)
    x , y = rand(1:n) , rand(1:n) # Generating the postion randomically
    pos = n*(x-1)+y # Converting (x,y) to the position on the Vertices Array

    # Calculating neighbors' spin sum
    spinclose = sum(spinneighbors(spinmatrix, x, y))

    # Calculating the energy on this environment
    energy_up = -spinclose - h
    energy_down = spinclose + h

    # Calculating the spin new "probabilities"
    prob_up = exp(-energy_up/temp)
    prob_down = exp(-energy_down/temp)
    normalization = prob_up + prob_down

    # Changing spin with the new probabilities
    towersort = rand()*normalization
    if towersort < prob_up
        spinmatrix[i][j] = 1
        # println("Changed spin ($(int(x)),$(int(y))) to 1")
    else
        spinmatrix[i][j] = -1
        # println("Changed spin ($(int(x)),$(int(y))) to -1")
    end
end

# HeatBath Algorithm for a Given Graph at a Given Temperature
function heatbath!(spinmatrix::Array{Int64,2}; temp::Float64=1.0, h::Float64=0.0, maxit::Int=20000, plot::Bool=true)
    xi = 1:maxit
    mi = Array(Float64, 0)
    
    for iter in xi
        heatbathstep!(spinmatrix, temp=temp,h=h)
        push!(mi, abs(magnetization(spinmatrix)))
    end
    
    println("Finished with magnetization $(mi[end])")
    if plot
        PyPlot.plot(xi, mi, "-", color="blue")
        PyPlot.ylim(0,1.2)
        PyPlot.title("HeatBath on Ising for T=$temp")
        PyPlot.xlabel("Number of Iterations")
        PyPlot.ylabel("Magnetization")
       
        PyPlot.savefig("Plots/HeatBath/hb_mag_$temp.png")
        PyPlot.close()
        # plot(x=xi,y=mi, Geom.point, Geom.line) #Gadfly
    end
    return mi[end]
end

# HeatBath Algorithm for many Graphs at a Given Temperature
function heattemp(n::Int, temp::Float64; qtd::Int=100, h::Float64=0.0, maxit::Int=20000)
    mi = Array(Float64,0)
    for i in 1:qtd
        g = spingrid(n)
        push!(mi, heatbath!(g,temp=temp,h=h,maxit=maxit, plot=false))
    end
    return sum(mi)/qtd # or just mean(mi) ? 
end

# HeatBath Algorithm for many Graphs at several Temperatures
function heatrange(n::Int; maxtemp::Float64=6.0, qtd::Int=100, h::Float64=0.0, maxit::Int=20000)
    mα = Array(Float64,0)
    for i in 0.1:0.1:maxtemp
       push!(mα,heattemp(n, i, qtd=qtd, h=h, maxit=maxit))
    end

    PyPlot.plot(0.1:0.1:maxtemp, mα, "-", color="red")
    PyPlot.title("Magnetization over Temperatures with HeatBath")
    PyPlot.savefig("Plots/heatbath_$(int(maxtemp))")
    PyPlot.close()
end