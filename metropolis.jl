# Step of the Metropolis' Algorithm
function metropolisstep!(g::GenericGraph; temp::Float64 = 1.0, h::Float64=0.0)
    if temp == 0 error("The Temperature can't be ZERO!!") end

    # Inicializations
    n = sqrt(num_vertices(g))
    x , y = rand(1:n) , rand(1:n) # Generating the postion randomically
    pos = n*(x-1)+y # Converting (x,y) to the position on the Vertices Array
    pos_spin = g.vertices[pos].attributes["spin"]

    # Calculating neighbors' spin sum
    mclose = sum([neighbor.attributes["spin"] for neighbor in out_neighbors(vertices(g)[pos],g)])

    # Calculating the energy on this environment
    energy_up = -mclose - h
    energy_down = mclose + h
    ΔE = (energy_up - energy_down)*(pos_spin) # energy necessary to flip the spin
    if ΔE < 0
        pos_spin *= -1 # flip spin
    else 
        # ELSE flip with Metropolis probability
        prob_flip = exp(-ΔE/temp)
        if rand() < prob_flip
            pos_spin *= -1 # flip spin
        end
    end
end

# Metropolis' Algorithm
function metropolis!(g::GenericGraph; temp::Float64=1.0, h::Float64=0.0, maxit::Int=10000, plot::Bool=true)
    xi = 1:maxit
    mi = Array(Float64, 0)
    
    for iter in xi
        metropolisstep!(g, temp=temp,h=h)
        push!(mi, abs(magnetization(g)))
    end
    
    println("Finished with magnetization $(mi[end])")
    if plot
	    PyPlot.plot(xi, mi, "o", color="blue")
	    PyPlot.plot(xi, mi, "-", color="blue")
	    PyPlot.ylim(0,1.2) # nao existia
	
	    PyPlot.title("Metropolis on Ising for T=$temp")
	    PyPlot.xlabel("Number of Iterations")
	    PyPlot.ylabel("Magnetization")
    	
	    #PyPlot.savefig("Plots/Metropolis/metro_mag_$temp.png")
	    #PyPlot.close()
	    # plot(x=xi,y=mi, Geom.point, Geom.line) #Gadfly
	end
    return mi[end]
end