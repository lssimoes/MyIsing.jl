##########################################
#                                        #
#   Wolff Functions Implementation       #
#                                        #
##########################################

# Step (RECURSIVE) of Wolff's Algorith
function wolffclusterstep!(spinmatrix::Array{Int, 2},
                           cluster::BitArray{2},
                           i::Int, 
                           j::Int; 
                           temp::Float64 = 1.0, 
                           h::Float64    = 0.0)

    if temp == 0 error("The Temperature can't be ZERO!!") end
    cluster[i,j] = true

    for (ni, nj) in neighbors(spinmatrix, i, j)
        if cluster[ni,nj] == false && spinmatrix[ni,nj] == spinmatrix[i,j] #dunno if the order matters
           if rand() < 1 - exp(-2/temp) wolffclusterstep!(spinmatrix, cluster, ni, nj, temp=temp, h=h) end
        end
    end
end

# Function that calculates the probabillity o fflipping the cluster given a external field h
wolffflipprob(grid::Array{Int, 2}, clt::BitArray{2}, h::Float64) = (prob=exp(-2h*clusterspin(grid, clt))) < 0 ? 1:prob

# Wolff's Algorithm for a Given Graph at a Given Temperature
function wolff!(spinmatrix::Array{Int, 2};
                temp::Float64 = 1.0,
                h::Float64    = 0.0,
                maxit::Int    = 300,
                plot::Bool    = true,
                verbose::Bool = true)

    # Defined here because I need this variable inside 'verbose' clause
    mi = Float64[]

    for iter in 1:maxit
        # Setting variables
        x , y = rand(1:size(spinmatrix, 1)) , rand(1:size(spinmatrix, 12))
        cluster = falses(size(spinmatrix))

        # Calculating Wolff step
        wolffclusterstep!(spinmatrix, cluster, x, y, temp=temp, h=h)
        if h == 0.0 || rand() < wolffflipprob(spinmatrix, cluster, h) flip!(spinmatrix, cluster) end

        # Appending the result
        push!(mi, abs(magnetization(spinmatrix)))
    end

    if verbose 
        # Saving to a .csv that informs Method, Size, Temperature and QtdIterations
        df = DataFrame(Iterations=1:maxit,Magnetization=mi)
        pathcsv = "Data/Wolff/wolff_$(size(spinmatrix,1))grid_$(temp)temp_$(int(h))h_$(maxit)iterations"
        writetable(pathcsv, df)

        println("Finished with magnetization $(mi[end])") 
        println("Data saved to file: " * pathcsv)
    end
    if plot plotcsv(pathcsv) end

    return mi[end]
end