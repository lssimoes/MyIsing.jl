# Step of the HeatBath Algorithm
function heatbathstep!(g::GenericGraph; temp::Float64 = 1.0, h::Float64=0.0)
    if temp == 0 error("The Temperature can't be ZERO!!") end

    # Inicializations
    n = sqrt(num_vertices(g))
    x , y = rand(1:n) , rand(1:n) # Generating the postion randomically
    pos = n*(x-1)+y # Converting (x,y) to the position on the Vertices Array

    # Calculating neighbors' spin sum
    spinclose = 0
    for neighbor in out_neighbors(vertices(g)[pos],g)
        spinclose += neighbor.attributes["spin"]
    end

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
        g.vertices[pos].attributes["spin"] = 1
        # println("Changed spin ($(int(x)),$(int(y))) to 1")
    else
        g.vertices[pos].attributes["spin"] = -1
        # println("Changed spin ($(int(x)),$(int(y))) to -1")
    end
end

# HeatBath Algorithm for a Given Graph at a Given Temperature
function heatbath!(g::GenericGraph; temp::Float64=1.0, h::Float64=0.0, maxit::Int=10000)
    xi = 1:maxit
    mi = Array(Float64, 0)
    
    for iter in xi
        heatbathstep!(g, temp=temp,h=h)
        push!(mi, abs(magnetization(g)))
    end
    
    println("Finished with magnetization $(mi[end])")
    #PyPlot.plot(xi, mi, "-", color="blue")
    #PyPlot.ylim(0,1.2) # nao existia

    #PyPlot.title("HeatBath on Ising for T=$temp")
    #PyPlot.xlabel("Number of Iterations")
    #PyPlot.ylabel("Magnetization")
    
    #PyPlot.savefig("Plots/HeatBath/hb_mag_$temp.png")
    #PyPlot.close()
    # plot(x=xi,y=mi, Geom.point, Geom.line) #Gadfly
    return mi[end]
end

# HeatBath Algorithm for many Graphs at a Given Temperature
function heattemp(n::Int, temp::Float64; qtd::Int=10, h::Float64=0.0, maxit::Int=50000)
    mi = Array(Float64,0)
    for i in 1:qtd
        g = nsquaregraph(n)
        push!(mi, heatbath!(g,temp=temp,h=h,maxit=maxit))
    end
    return sum(mi)/qtd
end

# HeatBath Algorithm for a set Graph at several Temperatures
# TODO HeatBath Algorithm for many Graphs at several Temperatures
function heatrange(n::Int; h::Float64=0.0, maxtemp::Float64=6.0, maxit::Int=10000 )
    ti = 0.1:0.1:maxtemp
    mi = Array(Float64, 0)
    arr = spinmatrix(n)
    
    for temp in ti
        g = matrix2graph(arr)
        push!(mi, heatbath!(g,temp=temp, h=h, maxit=maxit))
    end

    PyPlot.plot(ti,mi, "-", color="red")
    PyPlot.title("Magnetization over Temperatures with HeatBath")
    PyPlot.savefig("Plots/heatbath_$(int(maxtemp))")
    PyPlot.close()
end