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
randspin()                                = rand(0:1)*2-1
# To calculate the magnetization of the system
magnetization(spinmatrix::Array{Int,2})   = sum([spin for spin in spinmatrix])/size(spinmatrix)[1]^2
# To generate a random spinmatrix
spinmatrix(n::Int)                        = [randspin() for i in 1:n, j in 1:n]

#   FIX IT ALL !!!!
#
#   julia> spinneighbors(6,arr)
#   4-element Array{Int64,1}:
#    1
#    1
#    1 # Extra element !!
#    1
#   
#   julia> arr
#   3x3 Array{Int64,2}:
#     1  -1  1
#    -1   1  1
#     1  -1  1



# Creates a list of the POS neighbors' positions
function neighbors(pos::Int, n::Int)
    actual = Array(Int, 0)
    possible = [pos+1 , pos-n , pos-1 , pos+n]
    println(possible)
    for elem in possible
        if elem > 0 && elem > 0 && elem <= n^2 && elem <= n^2
            push!(actual, elem)
        end
    end
    println(actual)
    return actual
end
# Returns the spins of POS heighbors at spinmatrix
function spinneighbors(pos::Int, spinmatrix::Array{Int,2})
    n = size(spinmatrix)[1]
    spins = Array(Int, 0)
    for i in neighbors(pos, n)
        push!(spins, spinmatrix[i])
        println(spins)
        println("Passei")
    end
    return spins
end

#   EXPLAINING FUNCTIONALITY
#   TO "FIX" JULIA NUMBERING
#   OF MATRICES
#
#   julia> [1 2 3; 4 5 6;7 8 9]
#   3x3 Array{Int64,2}:
#    1  2  3
#    4  5  6
#    7  8  9
#   
#   julia> [1 2 3; 4 5 6;7 8 9].'
#   3x3 Array{Int64,2}:
#    1  4  7
#    2  5  8
#    3  6  9
#   
#   julia> [1 2 3; 4 5 6;7 8 9][4]
#   2
#   
#   julia> [1 2 3; 4 5 6;7 8 9].'[4]
#   4


end # Module END