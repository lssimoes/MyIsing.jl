##########################################
#                                        #
#   Plotting Functions Implementation    #
#                                        #
##########################################

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

    # TO FIX: Since wolff converges quicklier, we might changes the values a little bit
    if f == wolff!
        size *= 10
        ensembles /= 2
        maxit /= 100
    end

    for i in it
       push!(mα,meanontemp(f, size, i, qtd=ensembles, h=h, maxit=maxit))
    end

    if plot
        PyPlot.plot(it, mα, "-", color="red")
        PyPlot.title("Magnetization over Temperatures with " * "$f"[1:end-1])
        PyPlot.savefig("Plots/Phase/pyplot_" * "$f"[1:end-1] * "_$(size)grid_$(int(maxtemp))temp")
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