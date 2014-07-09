# Created by Lucas Silva Simões
# github.com/kaslusimoes
# MIT License

module Ising

using PyPlot, DataFrames, Gadfly


export randspin, 
       magnetization, 
       spingrid,
       flip!, 
       neighbors,
       spinneighbors,
       heatbath!,
       metropolis!,
       wolff!,
       phasediag,
       plotcsv

include("src/general.jl")
include("src/plotting.jl")
include("src/heatbath.jl")
include("src/metropolis.jl")
include("src/wolff.jl")

end