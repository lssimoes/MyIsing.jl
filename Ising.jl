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
       heatbath!,
       metropolis!,
       wolff!,
       phasediag


include("heatbath.jl")
include("metropolis.jl")
include("wolff.jl")

####################################################################
#                                                                  #
# Generic Functions that independ of the convergence algorithm     #
#                                                                  #
####################################################################

randspin()                                               = 2rand(0:1)-1                         # To sort random spins
magnetization(spinmatrix::Array{Int,2})                  = abs(mean(spinmatrix))                # To calculate the total magnetization
spingrid(n::Int)                                         = [randspin() for i in 1:n, j in 1:n]  # To generate a random spinmatrix
flip!(spinmatrix::Array{Int, 2}, i::Int, j::Int)         = spinmatrix[i, j] *= -1               # To flip a given spin
spinneighbors(spinmatrix::Array{Int, 2}, i::Int, j::Int) = [spinmatrix[i,j] for (i,j) in neighbors(spinmatrix, i, j)]

function neighbors(spinmatrix::Array{Int, 2}, i::Int, j::Int)
    n = Array((Int, Int), 0)

    if i > 1 push!(n, (i-1, j)) end
    if j > 1 push!(n, (i, j-1)) end

    if i < size(spinmatrix, 1) push!(n, (i+1, j)) end
    if j < size(spinmatrix, 2) push!(n, (i, j+1)) end

    return n
end

# Generic Algorithm that calculates a Function for many Graphs at several Temperature
function phasediag(f::Function; 
                   size::Int        = 10,
                   ensembles::Int   = 200,
                   h::Float64       = 0.0,
                   mintemp::Float64 = 0.2,
                   step::Float64    = 0.2,
                   maxtemp::Float64 = 6.0,
                   maxit::Int       = 50000,
                   plot::Bool       = true)

    mα = Array(Float64,0)
    it =  mintemp:step:maxtemp

    for i in it
       push!(mα,meanontemp(f, size, i, qtd=ensembles, h=h, maxit=maxit))
    end

    if plot
        PyPlot.plot(it, mα, "-", color="red")
        PyPlot.title("Magnetization over Temperatures with " * "$f"[1:end-1])
        PyPlot.savefig("Plots/" * "$f"[1:end-1] * "_$(size)grid_$(int(maxtemp))")
        PyPlot.close()
    end

    return it, mα
end

# Generic Algorithm that calculates a Function for many Graphs at a given Temperature
function meanontemp(f::Function, 
                    n::Int,
                    temp::Float64;
                    qtd::Int       = 200,
                    h::Float64     = 0.0,
                    maxit::Int     = 50000)

    mi = Array(Float64,0)

    for i in 1:qtd
        ensemble = spingrid(n)
        push!(mi, f(ensemble, temp=temp, h=h, maxit=maxit, plot=false, verbose=false))
    end
    println("Finished temperature $temp")

    return mean(mi) 
end

end # Module END