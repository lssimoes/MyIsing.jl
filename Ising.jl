# Created by Lucas Silva Simões
# github.com/kaslusimoes
# MIT License

module Ising

using Reexport, PyPlot #, Gadfly
@reexport using Graphs

export randspin, 
       magnetization, 
       spinmatrix,
       nsquaregraph, 
       matrix2graph,
       graph2matrix,
       hamiltonian, 
       heatbathstep!, 
       heatbath!,
       heatrange,
       metropolisstep!,
       metropolis!,
       wolffstep!

# To sort random spins
randspin()                     = rand(0:1)*2-1
# To calculate the magnetization of the system
magnetization(g::GenericGraph) = sum([vertex.attributes["spin"] for vertex in vertices(g)])/num_vertices(g)
# To generate a random spinmatrix
spinmatrix(n::Int)             = [randspin() for i in 1:n, j in 1:n]

# To create an n-square grid graph 
function nsquaregraph(n::Int)
    # Creating Vertices Array
    newvertices = Array(ExVertex, 0)
    for i in 1:n^2
        vertice = ExVertex(i,"eletron")
        vertice.attributes["spin"] = randspin()
        newvertices = [newvertices, vertice]
    end
    
    # Creating Edges Array
    n_edge = 4(n-2)^2 + 3*4(n-2) + 4*2
    n_edge = int(n_edge/2)
    newedges = Array(Edge{typeof(newvertices[1])}, 0)
    index = 1
    for i in 1:n-1
        for j in 1:n-1
            whoami = n*(i-1)+j
            newedges = [newedges, Edge(index, newvertices[whoami], newvertices[whoami+n])     # down neighbors
                                , Edge(index+1, newvertices[whoami], newvertices[whoami+1])]  # right
            index += 2
        end
        newedges = [newedges, Edge(index, newvertices[n*i], newvertices[n*i+n])] # down neighbor of the elemnt (i,n)
        index += 1
        # exception to treat last line, i represents the column instead of the line
        newedges = [newedges, Edge(index, newvertices[n*(n-1)+i], newvertices[n*(n-1)+i+1])] # neighbor to the right
        index += 1
    end
    return graph(newvertices, newedges, is_directed=false)
end

# To convert a Spin Matrix to a Spin Graph Grid
function matrix2graph(spinmatrix::Array{Int64,2})
    if size(spinmatrix,1) != size(spinmatrix,2) error("The SpinMatrix must be a Square matrix!") end
    n = size(spinmatrix)[1]
    newvertices = Array(ExVertex, 0)
    for i in 1:n
        for j in 1:n
            vertice = ExVertex(n*(i-1)+j, "eletron")
            vertice.attributes["spin"] = spinmatrix[i,j]
            newvertices = [newvertices, vertice]
        end
    end

    # Creating Edges Array
    n_edge = 4(n-2)^2 + 3*4(n-2) + 4*2
    n_edge = int(n_edge/2)
    newedges = Array(Edge{typeof(newvertices[1])}, 0)
    index = 1
    for i in 1:n-1
        for j in 1:n-1
            whoami = n*(i-1)+j
            newedges = [newedges, Edge(index, newvertices[whoami], newvertices[whoami+n])     # down neighbors
                                , Edge(index+1, newvertices[whoami], newvertices[whoami+1])]  # right
            index += 2
        end
        newedges = [newedges, Edge(index, newvertices[n*i], newvertices[n*i+n])] # down neighbor of the elemnt (i,n)
        index += 1
        # exception to treat last line, i represents the column instead of the line
        newedges = [newedges, Edge(index, newvertices[n*(n-1)+i], newvertices[n*(n-1)+i+1])] # neighbor to the right
        index += 1
    end
    return graph(newvertices, newedges, is_directed=false)
end

# (Paulo) To convert a Spin Graph Grid to a Spin Matrix
function graph2matrix(g::GenericGraph)
    n = int(sqrt(num_vertices(g)))
    spinarray = Array(Int, n, n)
    for i in 1:n
        for j in 1:n
            spinarray[i,j] = vertices(g)[n*(i-1) + j].attributes["spin"]
        end
    end

    return spinarray
end

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
function heatbath!(g::GenericGraph; temp::Float64=1.0, h::Float64=0.0, maxit::Int=10000, ε::Float64=0.001)
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
function heattemp(n::Int, temp::Float64; qtd::Int=10, h::Float64=0.0, maxit::Int=50000, ε::Float64=0.001)
    mi = Array(Float64,0)
    for i in 1:qtd
        g = nsquaregraph(n)
        push!(mi, heatbath!(g,temp=temp,h=h,maxit=maxit,ε=ε))
    end
    return sum(mi)/qtd
end

# HeatBath Algorithm for a set Graph at several Temperatures
# TODO HeatBath Algorithm for many Graphs at several Temperatures
function heatrange(n::Int; h::Float64=0.0, maxtemp::Float64=6.0, maxit::Int=10000, ε::Float64=0.001 )
    ti = 0.1:0.1:maxtemp
    mi = Array(Float64, 0)
    arr = spinmatrix(n)
    
    for temp in ti
        g = matrix2graph(arr)
        push!(mi, heatbath!(g,temp=temp, h=h, maxit=maxit, ε=ε))
    end

    PyPlot.plot(ti,mi, "-", color="red")
    PyPlot.title("Magnetization over Temperatures with HeatBath")
    PyPlot.savefig("Plots/heatbath_$(int(maxtemp))")
    PyPlot.close()
end

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
    else 
        # ELSE flip with Metropolis probability
        prob_flip = exp(-ΔE/temp)
        towersort = rand()
        if towersort < prob_flip
            g.vertices[pos].attributes["spin"] *= -1
        end
    end
end

# Metropolis' Algorithm
function metropolis!(g::GenericGraph; temp::Float64=1.0, h::Float64=0.0, maxit::Int=10000, ε::Float64=0.001)
    xi = 1:maxit
    mi = Array(Float64, 0)
    
    for iter in xi
        heatbathstep!(g, temp=temp,h=h)
        push!(mi, abs(magnetization(g)))
    end
    
    println("Finished with magnetization $(mi[end])")
    #PyPlot.plot(xi, mi, "o", color="blue")
    PyPlot.plot(xi, mi, "-", color="blue")
    PyPlot.ylim(0,1.2) # nao existia

    PyPlot.title("Metropolis on Ising for T=$temp")
    PyPlot.xlabel("Number of Iterations")
    PyPlot.ylabel("Magnetization")
    
    PyPlot.savefig("Plots/Metropolis/metro_mag_$temp.png")
    PyPlot.close()
    # plot(x=xi,y=mi, Geom.point, Geom.line) #Gadfly
    return mi[end]
end

# Step (RECURSIVE) fo Wolff's Algorith
function wolffstep!(pos::Int, g::GenericGraph; temp::Float64=1.0, h::Float64=0.0, maxit::Int=10000, ε::Float64=0.001)
    if temp == 0 error("The Temperature can't be ZERO!!") end

    pos_spin = g.vertices[pos].attributes["spin"]

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
            wolffstep!(pos,g,temp=temp, h=h, maxit=maxit, ε=ε)
        end
    end
end

# Wolff's Algorithm
function wolff!(g::GenericGraph; temp::Float64=1.0, h::Float64=0.0, maxit::Int=10000, ε::Float64=0.001)
    # Inicializations
    n = sqrt(num_vertices(g))
    x , y = rand(1:n) , rand(1:n) # Generating the postion randomically
    pos = n*(x-1)+y # Converting (x,y) to the position on the Vertices Array 

    wolffstep!(pos,g,temp=temp, h=h, maxit=maxit, ε=ε)
end

end # Module END