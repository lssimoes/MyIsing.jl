##########################################
#                                        #
#   Plotting Functions Implementation    #
#                                        #
##########################################

# Function that plots PyPlot and Gadfly given a path to a CSV
function plotcsv(path::ASCIIString)
    df = readtable(path)
    method, n, var, maxit = extract(path) # 'var' might be QtdEnsembles or Temperature

    if "$(names(df)[1])" == "Iterations" # then 'var' means Temperature
        label = "size $(n) temp $(var) maxit $(maxit)" 
    else #if "$(names(df)[1])" == "Temperature" #then 'var' means QtdEnsembles
        label = "size $(n) ensembles $(var) maxit $(maxit)" 
    end
    
    titleplot = "$(names(df)[1]) over $(names(df)[2]) using $(method),\n $(label)"

    plt = Gadfly.plot(df, x="$(names(df)[1])",y="$(names(df)[2])", Geom.point, Guide.xticks(ticks=[0:0.5:10.0]), Guide.title(titleplot))
    draw(PDF("Plots/$(method)/$(label).pdf",6inch,4inch),plt)
    draw(PNG("Plots/$(method)/$(label).png",6inch,4inch),plt)
    draw(SVGJS("Plots/$(method)/$(label).js.svg",6inch,4inch),plt)
end

# Function that returns the parameters given at a CSV path
# TEST: path = "Data/bla_23ble_17bli_42blu"
function extract(path::ASCIIString)
    pathargs = split(split(path, "/")[end], "_")

    method = "$(uppercase(pathargs[1][1]))$(pathargs[1][2:end])"
    n = match(r".[0-9]", "$(pathargs[2])").match
    arg2 = match(r".*[0-9]", "$(pathargs[3])").match
    maxit = match(r".*[0-9]", "$(pathargs[4])").match

    return method, n, arg2, maxit
end

# Generic Algorithm that calculates a Function for many Graphs at several Temperature
function phasediag(f::Function; 
                   n::Int        = 10,
                   ensembles::Int   = 200,
                   h::Float64       = 0.0,
                   mintemp::Float64 = 0.2,
                   step::Float64    = 0.2,
                   maxtemp::Float64 = 6.0,
                   maxit::Int       = 50000,
                   plot::Bool       = true)

    mi = Array(Float64,0)
    ti =  mintemp:step:maxtemp

    # TO FIX: Since wolff converges quicklier, we might changes the values a little bit
    if f == wolff!
        n *= 10
        ensembles = int(ensembles/4)
        maxit = int(maxit/200)
    end

    for i in ti
       push!(mi,meanontemp(f, n, i, qtd=ensembles, h=h, maxit=maxit))
    end

    # Saving to a .csv that informs Method, Size, QtdEnsembles and QtdIterations
    df = DataFrame(Temperature=ti,Magnetization=mi)
    method = "$(uppercase("$f"[1]))"*"$f"[2:end-1]

    pathcsv = "Data/$(method)/$(n)grid_$(ensembles)ensembles_$(maxit)iterations.csv"
    writetable(pathcsv, df)
    println("Data saved to file: " * pathcsv)

    if plot plotphase(pathcsv) end

    return df
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