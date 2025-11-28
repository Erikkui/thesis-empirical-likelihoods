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
# using BenchmarkTools: @btime

include("include_lib.jl")
include("wrapper.jl")
include("OU_solve.jl")
include("OU_main.jl")
include("create_synth_data.jl")


#####   Include functions in lib, set paths
CASE_NAME = "OU"
INCLUDE_PATH = joinpath( pwd(), "lib" )  # Create absolute path to lib
FIG_PATH = joinpath( pwd(), CASE_NAME, "figures" )  # Path to save figures
NETCDF_PATH = joinpath( pwd(), CASE_NAME, "netcdf") # Path to save NetCDF files
include_lib( INCLUDE_PATH )  # Recursively includes all functions in lib folder
#####


function run_OU()
    #### DATA LOADING AND SAVING OPTIONS ####
    synthetic_data = true       # Create and use synthetic data
    save_figures = false        # Save figures
    show_figures = true         # Show figures
    save_netcdf_file = true     # Save NetCDF file
    # datafile = "blowflies.csv"  # Data file name, used only if synthetic_data = false

    # Initial parameters for OU model, also used to generate synthetic data
    thetaa = 1.2
    sigma = 0.75

    theta = [ thetaa, sigma ]

    #### OPTIONS FOR RUN ####
    ## Which dataset to use for blowflies.csv
    dataset = 2

    ## Methods options
    use_diff = 1
    diff_order = [1]  # Orders of differences to calculate, e.g., [1, 2] for first and second order
    use_2D = false  # Use 2D array  [R0; R0_diff] instead of 1D array [R0, R0_diff] for R0 and R0_diff
    case = "bsl"  # "bsl" or "gsl"
    C_how = "cov"  # "cov" or "don" for standard covariance or Donsker theorem covariance
    axis_unif = "xax"  # "xax", "yax", or "log"
    use_log = "nolog"  # "log" or "nolog" for log transform for summary statistics

    ## Summary statistics calculation options
    eCDF = 1  # 0 for no eCDF, 1 for eCDF
    LL = [-2, 0]  # 0 for CIL, 1, 2... for ID; -1 to compute ecdf from signal too in addition to distances; -2 to not use distances
    chamfer = 0  # 0 for no chamfer distance, 1 for chamfer distance
    chamfer_k = [1, 2, 3] # Neighbors to consider for chamfer distance
    resample = 1 # 0 for no resampling, 1 for resampling (possibly to be deprecated)
    nsim = 30   # Number of model simulations per proposal theta (GSL: nsim = 1)
    nrep = 10  # Number of resamplings from simulations (always > 1)
    nbin = 10  # Number of bins for summary statistics

    ## Resampling options (BSL: bins; GSL: bins and data cov/mean)
    res_nrep = 25   # GSL only: res_nrep*res_nsamp iterations for data cov/mean calculation
    res_nsamp = 20  # Number of resamples for bin calc (BSL/GSL)

    #### MCMC OPTIONS ####
    nsimu = 10000   # MCMC chain length
    update_int = 30
    adapt_int = 20
    npar = length(theta)
    qcov = Matrix{Float64}( I, npar, npar )*1e-3

    sigma2 = 1.0
    noise = 0.0
    ssfun = like_eval

    #### SYNTHETIC DATA OPTIONS #### (only needed if synthetic_data = true)
    N_init = 10
    dt = 0.05
    N_end = 10
    nepo = 10

    # Save options in dictionaries
    data = Dict{Symbol, Any}(
        :synthetic_data => synthetic_data,
        :save_figures => save_figures,
        :show_figures => show_figures,
        :save_netcdf => save_netcdf_file,
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
        :res_dim => [res_nrep, res_nsamp],
        :params => theta,
        :synth_init => N_init,
        :synth_dt => dt,
        :nepo => nepo,
        :N_end => N_end,
        :minmax => nothing,  # Is filled later
        :nL => sum( LL .>= -1 ),
        :use_2D => use_2D,
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

    chain, sschain, results, data = OU_main( data, mcmc_options, mcmc_models );
    println( "Run completed. Acceptance rate: ", results[:accept] )

    return chain, sschain, results, data
end

chain, _, results, data = run_OU()
