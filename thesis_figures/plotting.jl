using Distributions
using Statistics
using StatsBase
using StaticArrays
using DifferentialEquations

include("plots_functions.jl")
include("blowfly_solve.jl")
include("OU_solve.jl")
include("lorenz3.jl")

# Uncomment the lines for the desired model to analyze

## OU params
# case = "OU"
# dir_path = "/home/eki/GitHub/thesis-empirical-likelihoods/OU/netcdf"
# nc_files = filter(f -> endswith(f, ".nc"), readdir(dir_path))
# nc_paths = joinpath.(dir_path, nc_files)
# param_names = [ "\\lambda", "\\sigma" ]
# params_true = [3.0, 2.0]
# burn_in = 10000
# n_samples = 1000


# L3 params
case = "L3"
dir_path = "/home/eki/GitHub/thesis-empirical-likelihoods/L3/netcdf"
nc_files = filter(f -> endswith(f, ".nc"), readdir(dir_path))
nc_paths = joinpath.(dir_path, nc_files)
param_names = [ "\\sigma", "\\rho", "\\beta" ]
params_true = [10.0, 28.0, 8/3]
burn_in = 10000
n_samples = 1


## Blowfly params
# case = "blowfly"
# dir_path = "/home/eki/GitHub/thesis-empirical-likelihoods/blowfly/netcdf"
# nc_files = filter(f -> endswith(f, ".nc"), readdir(dir_path))
# nc_paths = joinpath.(dir_path, nc_files)
# println( nc_paths )
# param_names = [ "\\delta", "P", "N_0", "\\sigma^2_p", "\\tau", "\\sigma^2_d" ]
# delta = 0.16
# P = 6.5
# N_0 = 400
# sigma_p = 0.1
# tau = 14
# sigma_d = 0.1
# params_true = [ delta, P, N_0, sigma_p, tau, sigma_d ]
# burn_in = 15000
# n_samples = 1000
##


for file in nc_paths  # Select specific files to process
    printname = split( file, "/" )[end]
    println( "Processing file: "*printname )

    number = split( file, "_" )[1]
    number = split( number, "/")[end]

    name = split( file, "_" )[ end ]
    name = split( name, ".")[1]
    name = number*"_"*case*"_result_"*name
    name = joinpath( pwd(), "thesis_figures", case, name )

    fig, ax = plot_normality_checks(file)
    display( fig )
    # save( name*"_normality.pdf", fig )

    fig, ax = plot_mcmc_results(file, :histchain, param_names, params_true; burn_in=burn_in)
    display( fig )
    # save( name*"_mcmc.pdf", fig )

    if case != "L3"
        fig = plot_model_predictions(file, case; burn_in=burn_in, n_samples=n_samples)
        display( fig )
        # save( name*"_predictions.pdf", fig )
    end

end

##
labels = string.(collect(1:length(nc_paths)))
fig = plot_forest_multi(nc_paths, param_names, params_true; labels = labels, burn_in=burn_in)
display( fig )
##
name = case*"_forest.pdf"
# save( joinpath( pwd(), "thesis_figures", case, name ), fig )
