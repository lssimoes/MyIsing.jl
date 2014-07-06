# Step of the Metropolis' Algorithm
function metropolisstep!(g::GenericGraph; temp::Float64 = 1.0, h::Float64=0.0)
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
    ΔE = (energy_up - energy_down)*(g.vertices[pos].attributes["spin"]) # energy necessary to flip the spin
    if (ΔE < 0)
        g.vertices[pos].attributes["spin"] *= -1 # flip spin
        println("\nflipei")
    else 
    	print("entrei")
        # ELSE flip with Metropolis probability
        prob_flip = (energy_up/energy_down)^(g.vertices[pos].attributes["spin"])
        towersort = rand()
        if towersort < prob_flip
        	println(" e flipei!")
            g.vertices[pos].attributes["spin"] *= -1
        else
        	println(" mas não fiplei..")
        end
    end
end

# Metropolis' Algorithm
function metropolis!(g::GenericGraph; temp::Float64=1.0, h::Float64=0.0, maxit::Int=10000)
    xi = 1:maxit
    mi = Array(Float64, 0)
    
    for iter in xi
        metropolisstep!(g, temp=temp,h=h)
        push!(mi, abs(magnetization(g)))
    end
    
    println("Finished with magnetization $(mi[end])")
    PyPlot.plot(xi, mi, "o", color="blue")
    PyPlot.plot(xi, mi, "-", color="blue")
    PyPlot.ylim(0,1.2) # nao existia

    PyPlot.title("Metropolis on Ising for T=$temp")
    PyPlot.xlabel("Number of Iterations")
    PyPlot.ylabel("Magnetization")
    
    #PyPlot.savefig("Plots/Metropolis/metro_mag_$temp.png")
    #PyPlot.close()
    # plot(x=xi,y=mi, Geom.point, Geom.line) #Gadfly
    return mi[end]
end