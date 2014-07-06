# Step (RECURSIVE) of Wolff's Algorith
function wolffstep!(pos::Int, g::GenericGraph; temp::Float64=1.0, h::Float64=0.0, maxit::Int=10000)
    if temp == 0 error("The Temperature can't be ZERO!!") end

    pos_spin = g.vertices[pos].attributes["spin"]

    # Calculating neighbors' spin sum
    mclose = sum([neighbor.attributes["spin"] for neighbor in out_neighbors(vertices(g)[pos],g)])

    ΔE = 2(h+mclose)*pos_spin # energy necessary to flip the spin
    if (ΔE < 0)
        g.vertices[pos].attributes["spin"] *= -1 # flip spin
    else 
        # ELSE flip with Metropolis probability
        prob_flip = exp(-ΔE/temp)
        towersort = rand()
        if towersort < prob_flip
            g.vertices[pos].attributes["spin"] *= -1
        end
    end

    # Recursive flip neighbor with same spin
    for neighbor in out_neighbors(vertices(g)[pos],g)
        if neighbor.attributes["spin"] == pos_spin
            wolffstep!(pos,g,temp=temp, h=h, maxit=maxit)
        end
    end
end

# Wolff's Algorithm
function wolff!(g::GenericGraph; temp::Float64=1.0, h::Float64=0.0, maxit::Int=10000, plot::Bool=true)
    # Inicializations
    n = sqrt(num_vertices(g))
    x , y = rand(1:n) , rand(1:n) # Generating the postion randomically
    pos = n*(x-1)+y # Converting (x,y) to the position on the Vertices Array 

    wolffstep!(pos,g,temp=temp, h=h, maxit=maxit)

    println("Finished with magnetization $(mi[end])")
    if plot
        #PyPlot.plot(xi, mi, "o", color="blue")
        PyPlot.plot(xi, mi, "-", color="blue")
        PyPlot.ylim(0,1.2) # nao existia
    
        PyPlot.title("Wolff on Ising for T=$temp")
        PyPlot.xlabel("Number of Iterations")
        PyPlot.ylabel("Magnetization")
        
        PyPlot.savefig("Plots/Wolff/wolff_mag_$temp.png")
        PyPlot.close()
        # plot(x=xi,y=mi, Geom.point, Geom.line) #Gadfly
    end
    return mi[end]
end
