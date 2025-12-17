using Distributions

include("plots_functions.jl")
include("blowfly_solve.jl")
include("OU_solve.jl")
include("lorenz3.jl")


## OU params
# case = "OU"
# dir_path = "/home/eki/GitHub/thesis-empirical-likelihoods/OU/netcdf"
# nc_files = filter(f -> endswith(f, ".nc"), readdir(dir_path))
# nc_paths = joinpath.(dir_path, nc_files)
# param_names = [ "\\lambda", "\\sigma" ]


## L3 params
# case = "L3"
# dir_path = "/home/eki/GitHub/thesis-empirical-likelihoods/L3/netcdf"
# nc_files = filter(f -> endswith(f, ".nc"), readdir(dir_path))
# nc_paths = joinpath.(dir_path, nc_files)
# param_names = [ "\\sigma", "\\rho", "\\beta" ]


# ## Blowfly
case = "blowfly"
dir_path = "/home/eki/GitHub/thesis-empirical-likelihoods/blowfly/netcdf"
nc_files = filter(f -> endswith(f, ".nc"), readdir(dir_path))
nc_paths = joinpath.(dir_path, nc_files)
println( nc_paths )
param_names = [ "P", "N_0", "\\delta", "\\sigma^2_p", "\\tau", "\\sigma^2_d" ]
burn_in = 15000
##


for file in nc_paths[ [2, 4, 5, 6, 7]]
    prinname = split( file, "/" )[end]
    println( "Processing file: "*prinname )
    name = split( file, "_" )[end]
    name = split( name, ".")[1]
    name = case*"_result_"*name
    name = joinpath( pwd(), "thesis_figures", case, name )

    fig, ax = plot_normality_checks(file)
    display( fig )
    save( name*"_normality.pdf", fig )

    fig, ax = plot_mcmc_results(file, :histchain, param_names, burn_in=burn_in)
    display( fig )
    save( name*"_mcmc.pdf", fig )

    fig, ax = plot_model_predictions(file, case; burn_in=burn_in, n_samples=100)
    display( fig )
    save( name*"_predictions.pdf", fig )

end
