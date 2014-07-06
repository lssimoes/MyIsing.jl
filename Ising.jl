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
#       metrotemp,
       wolffstep!


include("heatbath.jl")
include("metropolis.jl")
include("wolff.jl")

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

end # Module END