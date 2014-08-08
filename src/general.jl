##########################################
#                                        #
#   Generic Functions Implementation     #
#                                        #
##########################################

randspin()                                               = 2rand(0:1)-1                         # To sort random spins
magnetization(spinmatrix::Array{Int,2})                  = abs(mean(spinmatrix))                # To calculate the total magnetization
spingrid(n::Int)                                         = [randspin() for i in 1:n, j in 1:n]  # To generate a random spinmatrix
flip!(spinmatrix::Array{Int, 2}, i::Int, j::Int)         = spinmatrix[i, j] *= -1               # To flip a given spin
spinneighbors(spinmatrix::Array{Int, 2}, i::Int, j::Int) = [spinmatrix[i,j] for (i,j) in neighbors(spinmatrix, i, j)] # To sum the spin nearby

# Function that flips the sites (x,y) that are 'true' on 'cluster'
function flip!(spinmatrix::Array{Int, 2}, cluster::BitArray{2})
    for y in 1:size(cluster,2) for x in 1:size(cluster,1) 
        if cluster[x,y] flip!(spinmatrix,x,y) end
    end end
end

# Function that cleans the cluster
# function cleancluster(cluster::BitArray{2})
#     for y in 1:size(cluster,2) for x in 1:size(cluster,1) 
#         if cluster[x,y] cluster[x,y] = false end
#     end end
# end

# Function that returns the 'amount' of spin on a cluster
function clusterspin(spinmatrix::Array{Int, 2}, cluster::BitArray{2})
    spin = 0
    for y in 1:size(cluster,2) for x in 1:size(cluster,1) 
        if cluster[x,y] spin += spinmatrix[x,y] end
    end end
    return spin
end

# Function that evaluates some position neighbors and returns a list of their positions (tuple)
function neighbors(spinmatrix::Array{Int, 2}, i::Int, j::Int)
    n = (Int,Int)[] # same as Array((Int, Int), 0)

    if i > 1 push!(n, (i-1, j)) end
    if j > 1 push!(n, (i, j-1)) end

    if i < size(spinmatrix, 1) push!(n, (i+1, j)) end
    if j < size(spinmatrix, 2) push!(n, (i, j+1)) end

    return n
end

