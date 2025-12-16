using Statistics
using LinearAlgebra
using Distributions
using Random: randperm
using Distances: pairwise, Euclidean
using StatsBase
using LoopVectorization
using NearestNeighbors
using CairoMakie
using NCDatasets
using ProgressBars
using CSV
# using BenchmarkTools: @btime

include("include_lib.jl")
include("blowfly_solve.jl")
include("blowfly_main.jl")
include("create_synth_data.jl")
include("wrapper.jl")


#####   Include functions in lib, set paths
CASE_NAME = "blowfly"
INCLUDE_PATH = joinpath( pwd(), "lib" )  # Create absolute path to lib
FIG_PATH = joinpath( pwd(), CASE_NAME, "figures" )  # Path to save figures
NETCDF_PATH = joinpath( pwd(), CASE_NAME, "netcdf") # Path to save NetCDF files
include_lib( INCLUDE_PATH )  # Recursively includes all functions in lib folder
#####

function run_blowfly()
    #### DATA LOADING AND SAVING OPTIONS ####
    synthetic_data = true       # Create and use synthetic data
    save_figures = false        # Save figures
    show_figures = true         # Show figures
    save_netcdf_file = true     # Save NetCDF file
    datafile = "blowflies.csv"  # Data file name, used only if synthetic_data = false

    # Initial parameters for blowfly model, also used to generate synthetic data
    # These parameters are based on the paper https://doi.org/10.48550/arXiv.1411.4564
    delta = 0.16
    P = 6.5
    N_0 = 400
    sigma_p = 0.1
    tau = 14
    sigma_d = 0.1
    theta = [ delta, P, N_0, sigma_p, tau, sigma_d ]

    #### OPTIONS FOR RUN ####
    ## Which dataset to use for blowflies.csv
    dataset = 4  # 1-4, for data in blowflies.csv

    ## Methods options
    use_diff = 1
    diff_order = [ 1, 2 ]  # Orders of differences to calculate, e.g., [1, 2] for first and second order
    case = "bsl"  # "bsl" or "gsl"
    C_how = "cov"  # "cov" or "don" for standard covariance or Donsker theorem covariance
    axis_unif = "yax"  # "xax", "yax", or "log"
    use_log = "nolog"  # "log" or "nolog" for log transform for summary statistics

    ## Summary statistics calculation options
    eCDF = 1  # 0 for no eCDF, 1 for eCDF
    LL = [ -1 ]  # CIL: 0 for distances, -1 for signal; ID: positive integers for kNN distances
    chamfer = 0 # 0 for no chamfer distance, 1 for chamfer distance
    chamfer_k = [ 1, 2, 3, 8, 9, 10, 11, 18, 19, 20 ] # Neighbors to consider for chamfer distance
    nsim = 20  # Number of model simulations per proposal theta (GSL: usually nsim = 1)
    nrep = 100 # Number of resamplings from simulations (always > 1)
    nbin = 10  # Number of bins for summary statistics

    ## Resampling options (BSL: bins; GSL: bins and data cov/mean)
    resample = 1    # 0 for no resampling, 1 for resampling
    res_like = 300   # Iterations for data cov/mean calculation
    res_bins = 100  # Number of resamples for bin calc (BSL/GSL)
    window = 80

    #### MCMC OPTIONS ####
    nsimu = 10000   # MCMC chain length
    update_int = 15
    adapt_int = 20
    npar = length(theta)
    qcov = Matrix{Float64}( I, npar, npar )*1e-2

    sigma2 = 1.0
    noise = 0.0
    ssfun = like_eval

    #### SYNTHETIC DATA OPTIONS #### (only needed if synthetic_data = true)
    init = 180  # Initial condition for synthetic data
    N = 200     # Length of synthetic data
    burn_in = 0    # Time after which the data is considered "stable"
    t = N+1  # Time vector from 0 to N
    nepo = 1

    # Save options in dictionaries
    data = Dict{Symbol, Any}(
        :synthetic_data => synthetic_data,
        :save_figures => save_figures,
        :show_figures => show_figures,
        :save_netcdf => save_netcdf_file,
        :params => theta,
        :dataset_num => dataset,
        :use_diff => use_diff,
        :diff_order => diff_order,
        :case => case,
        :C_how => C_how,
        :uni => axis_unif,
        :log => use_log,
        :eCDF => eCDF,
        :LL => LL,
        :chamfer => chamfer,
        :chamfer_k => chamfer_k,
        :resample => resample,
        :case_dim => [nsim, nrep],
        :nbin => nbin,
        :res_dim => [res_like, res_bins],
        :synth_dt => 1.0,
        :synth_init => init,
        :synth_N => t,
        :synth_burn_in => burn_in,
        :synth_t => t,
        :nepo => nepo,
        :minmax => nothing,  # Is filled later
        :nL => length(LL),
        :window => window,
    )

    mcmc_options = Dict{Symbol, Any}(
        :nsimu => nsimu,
        :update_int => update_int,
        :adapt_int => adapt_int,
        :qcov => qcov,
    )

    mcmc_models = Dict{Symbol, Any}(
        :sigma2 => sigma2,
        :noise => noise,
        :ssfun => ssfun
    )

    chain, sschain, results, data = blowfly_main( data, mcmc_options, mcmc_models, datafile=datafile );
    # @profview blowfly_main( data, mcmc_options, mcmc_models, datafile=datafile );
    println( "Run completed. Acceptance rate: ", results[:accept] )

    return chain, sschain, results, data
end

chain, _, _, data = run_blowfly()
