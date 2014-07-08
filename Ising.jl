# Created by Lucas Silva Simões
# github.com/kaslusimoes
# MIT License

module Ising

using PyPlot #, Gadfly

export randspin, 
       magnetization, 
       spinmatrix,
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
randspin()                                         = rand(0:1)*2-1
# To calculate the magnetization of the system
magnetization(spinmatrix::Array{Int,2})            = mean([spin for spin in spinmatrix])
# To generate a random spinmatrix
spingrid(n::Int)                                   = [randspin() for i in 1:n, j in 1:n]
# To flip the spin at a given position
flip!(spinmatrix::Array{Int, 2}, i::Int, j::Int) = spinmatrix[i, j] *= -1

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