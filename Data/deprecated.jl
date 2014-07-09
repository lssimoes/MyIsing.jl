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

function pygraphs()
    # PyPlot
    PyPlot.plot(ti, mi, "-", color="red")
    PyPlot.title("Magnetization over Temperatures with " * "$f"[1:end-1])
    PyPlot.xlabel("Temperature")
    PyPlot.ylabel("Magnetization") 
    PyPlot.ylim(0,1.1)
    PyPlot.savefig("Plots/Phase/pyplot_" * "$f"[1:end-1] * "_$(n)grid_$(int(maxtemp))temp")
    PyPlot.close()

    PyPlot.plot(xi, mi, "o", color="blue")
    PyPlot.title("HeatBath on Ising for n=$(size(spinmatrix,1)) and T=$temp")
    PyPlot.xlabel("Number of Iterations")
    PyPlot.ylabel("Magnetization")   
    PyPlot.ylim(0,1.1)
    PyPlot.savefig("Plots/HeatBath/hb_n$(size(spinmatrix,1))_temp$temp.png")
    PyPlot.close()
end