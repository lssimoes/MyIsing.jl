##########################################
#                                        #
#   HeatBath Functions Implementation    #
#                                        #
##########################################

# Step of the HeatBath Algorithm
function heatbathstep!(spinmatrix::Array{Int64,2}; 
                       temp::Float64  = 1.0,
                       h::Float64     = 0.0,
                       verbose::Bool  = true)

    if temp == 0 error("The Temperature can't be ZERO!!") end

    n = size(spinmatrix, 1)
    x , y = rand(1:n) , rand(1:n) # Generating the position randomically

    # Calculating neighbors' spin sum
    magnear = sum(spinneighbors(spinmatrix, x, y))

    # Calculating the energy on this environment
    energy_up = -magnear - h
    energy_down = magnear + h

    # Calculating the spin new "probabilities"
    prob_up = exp(-energy_up/temp)
    normalization = prob_up + exp(-energy_down/temp)

    # Changing spin with the new probabilities
    spinmatrix[x,y] = rand()*normalization < prob_up ? 1 : -1
    if verbose 
        println("Changed spin ($(int(x)),$(int(y))) to $(spinmatrix[i,j])") 
    end
end

# HeatBath Algorithm for a Given Graph at a Given Temperature
function heatbath!(spinmatrix::Array{Int64,2};
                   temp::Float64  = 1.0,
                   h::Float64     = 0.0,
                   maxit::Int     = 20000,
                   plot::Bool     = true,
                   verbose::Bool  = true)

    xi = 1:maxit
    mi = Array(Float64, 0)
    
    for iter in xi
        heatbathstep!(spinmatrix, temp=temp, h=h, verbose=false)
        push!(mi, abs(magnetization(spinmatrix)))
    end
    
    if verbose 
        # Saving to a .csv that informs Method, Size, Temperature and QtdIterations
        df = DataFrame(Iterations=xi,Magnetization=mi)
        pathcsv = "Data/Heatbath/heatbath_$(size(spinmatrix,1))grid_$(temp)temp_$(int(h))h_$(maxit)iterations"
        writetable(pathcsv, df)

        println("Finished with magnetization $(mi[end])") 
        println("Data saved to file: " * pathcsv)
    end
    if plot plotcsv(pathcsv) end

    return mi[end]
end