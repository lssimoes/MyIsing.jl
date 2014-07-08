# Created by Lucas Silva Simões
# github.com/kaslusimoes
# MIT License

module Ising

using PyPlot #, Gadfly

export randspin, 
       magnetization, 
       spingrid,
       flip!, 
       spinneighbors,
       heatbathstep!, 
       heatbath!,
       heatrange,


include("heatbath.jl")
#include("metropolis.jl")
#include("wolff.jl")

randspin()                                       = 2rand(0:1)-1                         # To sort random spins
magnetization(spinmatrix::Array{Int,2})          = abs(mean(spinmatrix))                # To calculate the total magnetization
spingrid(n::Int)                                 = [randspin() for i in 1:n, j in 1:n]  # To generate a random spinmatrix
flip!(spinmatrix::Array{Int, 2}, i::Int, j::Int) = spinmatrix[i, j] *= -1               # To flip a given spin

# Returns the spins of the neighbors of (i,j) at spinmatrix
function spinneighbors(spinmatrix::Array{Int, 2}, i::Int, j::Int)
    n = Array(Int, 0)
    if i > 1 push!(n, spinmatrix[i-1, j]) end
    if j > 1 push!(n, spinmatrix[i, j-1]) end

    if i < size(spinmatrix, 1) push!(n, spinmatrix[i+1, j]) end
    if j < size(spinmatrix, 2) push!(n, spinmatrix[i, j+1]) end

    return n
end

end # Module END