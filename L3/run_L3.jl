using Statistics
using LinearAlgebra
using Distributions
using Random: randperm
using Distances: pairwise, Euclidean
using StatsBase
using DifferentialEquations
using LoopVectorization
using NearestNeighbors
using CairoMakie
using NCDatasets
using ProgressBars
# using BenchmarkTools: @btime

include("include_lib.jl")
include("wrapper.jl")
include("lorenz3.jl")
include("L3_main.jl")
include("create_synth_data.jl")


#####   Include functions in lib, set paths
CASE_NAME = "L3"
INCLUDE_PATH = joinpath( pwd(), "lib" )  # Create absolute path to lib
FIG_PATH = joinpath( pwd(), CASE_NAME, "figures" )  # Path to save figures
NETCDF_PATH = joinpath( pwd(), CASE_NAME, "netcdf") # Path to save NetCDF files
include_lib( INCLUDE_PATH )  # Recursively includes all functions in lib folder
#####


function run_L3()
    #### DATA LOADING AND SAVING OPTIONS ####
    save_figures = false        # Save figures
    show_figures = true         # Show figures
    save_netcdf_file = false     # Save NetCDF file

    # Initial parameters for L3 model, also used to generate synthetic data
    theta = [10.0, 28.0, 8/3]

    #### OPTIONS FOR RUN ####

    ## Methods options
    use_diff = 1
    diff_order = [1]  # Orders of differences to calculate, e.g., [1, 2] for first and second order
    case = "bsl"  # "bsl" or "gsl"
    C_how = "cov"  # "cov" or "don" for standard covariance or Donsker theorem covariance
    axis_unif = "xax"  # "xax", "yax", or "log"
    use_log = "log"  # "log" or "nolog" for log transform for summary statistics

    ## Summary statistics calculation options
    eCDF = 1  # 0 for no eCDF, 1 for eCDF
    LL = [ 0, 1, 2]  # CIL: 0 for distances, -1 for signal; ID: positive integers for kNN distances
    chamfer = 1  # 0 for no chamfer distance, 1 for chamfer distance
    chamfer_k = [1, 2] # Neighbors to consider for chamfer distance
    nsim = 10   # Number of model simulations per proposal theta (GSL: nsim = 1)
    nrep = 1  # Number of resamplings from simulations (always > 1)
    nbin = 10  # Number of bins for summary statistics

    ## Resampling options (BSL: bins; GSL: bins and data cov/mean)
    resample = 1    # 0 for no resampling, 1 for resampling for bins
    res_nrep = 100   # Iterations for data cov/mean calculation
    res_nsamp = 20  # Number of resamples for bin calc (BSL/GSL)

    #### MCMC OPTIONS ####
    nsimu = 1000   # MCMC chain length
    update_int = 30
    adapt_int = 50
    npar = length(theta)
    qcov = Matrix{Float64}( I, npar, npar )*1e-2

    sigma2 = 1.0
    noise = 0.0
    ssfun = like_eval

    #### SYNTHETIC DATA OPTIONS #### (only needed if synthetic_data = true)
    init = [12.577, 19.471, 23.073]
    dt = 1
    N_end = 300+10
    nepo = 1

    # Save options in dictionaries
    data = Dict{Symbol, Any}(
        :synthetic_data => true,
        :save_figures => save_figures,
        :show_figures => show_figures,
        :save_netcdf => save_netcdf_file,
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
        :res_dim => [res_nrep, res_nsamp],
        :params => theta,
        :synth_init => init,
        :synth_dt => dt,
        :nepo => nepo,
        :N_end => N_end,
        :minmax => nothing,  # Is filled later
        :nL => length(LL),
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

    # chain, sschain, results, data = L3_main( data, mcmc_options, mcmc_models );
    @profview L3_main( data, mcmc_options, mcmc_models )
    println( "Run completed. Acceptance rate: ", results[:accept], "\n" )

    return chain, sschain, results, data
end

chain, _, results, data = run_L3()
