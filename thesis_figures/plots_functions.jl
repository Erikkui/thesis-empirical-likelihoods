using CairoMakie
using NCDatasets
using LaTeXStrings
using MathTeXEngine
using Statistics


function plot_mcmc_results(nc_path::String, plot_type::Symbol, param_names::Vector{String}; burn_in::Int=0)
    """
    plot_mcmc_results(nc_path::String, plot_type::Symbol, param_names::Vector{String}; burn_in::Int=0)

    Plot MCMC results from NetCDF file.
    Now supports `:histchain` to plot chains and rotated histograms side-by-side.
    """

    # 1. Load Data
    ds = NCDataset(nc_path, "r")
    raw_chain = ds["chain"][:, :]

    # Read True Parameters (Robust Parsing)
    # if haskey(ds.attrib, "params")
    #     param_str = ds.attrib["params"]
    #     # Remove brackets [], commas, semicolons, then split
    #     clean_str = replace(param_str, r"[\[\],;]" => " ")
    #     true_params = parse.(Float64, split(clean_str))
    #     println("Parsed true parameters: ", true_params)
    # else
    #     error("Global attribute 'params' not found in NetCDF file.")
    # end
    true_params = raw_chain[1, :]
    close(ds)

    # 2. Apply Burn-in
    if burn_in >= size(raw_chain, 1)
        error("Burn-in length ($burn_in) exceeds chain length.")
    end
    chain = raw_chain[burn_in+1:end, :]
    n_samples, n_params = size(chain)

    if n_params != length(param_names)
        error("Parameter names count ($(length(param_names))) != data columns ($n_params).")
    end
    if length(true_params) != n_params
        error("Parsed true params ($(length(true_params))) != data columns ($n_params).")
    end

    # 3. Setup Theme
    colors = Makie.to_colormap(:seaborn_dark)
    set_theme!(palette = (; color = colors), fontsize = 24)

    # 4. Initialize Figure & Layout
    if plot_type == :chains
        n_cols = 2
        n_rows = ceil(Int, n_params / n_cols)
        fig_size = (900, 150 * n_rows)
        fig = Figure(size = fig_size, figure_padding = (3, 15, 3, 3))
    elseif plot_type == :histogram
        n_cols = ceil(Int, sqrt(n_params))
        n_rows = ceil(Int, n_params / n_cols)
        fig_size = (900, 300 * n_rows)
        fig = Figure(size = fig_size, figure_padding = (3, 15, 3, 3))
    elseif plot_type == :histchain
        # New Layout: 1 row per parameter, wider figure for side-by-side
        n_rows = n_params
        fig_size = (1100, 180 * n_rows)
        fig = Figure(size = fig_size, figure_padding = (5, 20, 5, 5))
    else
        error("Invalid plot_type. Options: :chains, :histogram, :histchain")
    end

    axs = Axis[]

    # 5. Plotting Loop
    for i in 1:n_params
        data_vec = chain[:, i]
        param_tex = latexstring("\\mathrm{" * param_names[i] * "}")

        # --- NEW: HistChain Logic ---
        if plot_type == :histchain
            # 1. Chain Axis (Left)
            ax_chain = Axis(fig[i, 1],
                ylabel = param_tex,
                # xlabel = (i == n_params) ? L"\mathrm{Iteration}" : "",
                xticklabelsvisible = (i == n_params),
                xgridvisible = true, ygridvisible = true,
                xlabelsize = 28, ylabelsize = 28,
                # yticks = Makie.WilkinsonTicks(4), # Same tick logic as chains
                # xticklabelsize = 18, yticklabelsize = 18
            )

            # 2. Histogram Axis (Right)
            ax_hist = Axis(fig[i, 2],
                ylabel = "",
                xticklabelsvisible = false, xticksvisible = false,
                yticklabelsvisible = false, yticksvisible = false,
                xgridvisible = true, ygridvisible = true
            )

            # Link Y axes so they share the parameter scale
            linkyaxes!(ax_chain, ax_hist)

            # Plot Data
            iterations = (burn_in+1):(burn_in+n_samples)
            lines!(ax_chain, iterations, data_vec, linewidth=1.5)

            # Rotated Density (direction=:x makes it horizontal)
            hist!(ax_hist, data_vec, normalization = :pdf, direction = :x, strokewidth = 1, color = colors[3], strokecolor = :black, bins = 20)

            # True Parameter Lines (Horizontal)
            hlines!(ax_chain, [true_params[i]], color = (colors[2], 0.85), linewidth = 3.75)
            hlines!(ax_hist, [true_params[i]], color = (colors[2], 0.85), linewidth = 3.75)

            # Set Column Ratios (80% chain, 20% histogram)
            colsize!(fig.layout, 1, Relative(0.8))
            colsize!(fig.layout, 2, Relative(0.2))

            push!(axs, ax_chain)
            push!(axs, ax_hist)

        else
            # --- EXISTING Logic for :chains and :histogram ---
            row = (plot_type == :chains) ? div(i-1, 2) + 1 : div(i-1, n_cols) + 1
            col = (plot_type == :chains) ? rem(i-1, 2) + 1 : rem(i-1, n_cols) + 1

            is_bottom = (plot_type == :chains) ? (i + 2 > n_params) : (i + n_cols > n_params)
            lbl_iter = L"\mathrm{Iteration}"
            lbl_dens = L"\mathrm{Density}"

            ticks_y = (plot_type == :chains) ? Makie.WilkinsonTicks(4) : Makie.automatic

            ax = Axis(fig[row, col],
                xlabel = (plot_type == :chains && is_bottom) ? lbl_iter : "",
                ylabel = (plot_type == :chains) ? param_tex : lbl_dens,
                # title = (plot_type == :histogram) ? param_tex : "",
                xgridvisible = true,
                ygridvisible = true,
                yticks = ticks_y,
                # xticklabelsize = 18,
                # yticklabelsize = 18
            )
            push!(axs, ax)

            if plot_type == :chains
                iterations = (burn_in+1):(burn_in+n_samples)
                lines!(ax, iterations, data_vec, linewidth=1.5)
                if !is_bottom
                    hidexdecorations!(ax, grid = false)
                end

            elseif plot_type == :histogram
                hist!(ax, data_vec, normalization = :pdf, color = (:black, 0.1), strokewidth = 1, strokecolor = :black)
                # density!(ax, data_vec, color = :transparent, strokewidth = 3.75, strokecolor = :black)
                vlines!(ax, [true_params[i]], color = (:red, 0.7), linewidth = 3.75)
            end
        end
    end

    # 6. Adjust spacing
    if plot_type == :histchain
        colgap!(fig.layout, 5) # Tight gap between chain and hist
        rowgap!(fig.layout, 10)
    else
        rowgap!(fig.layout, 10)
        colgap!(fig.layout, 20)
    end

    return fig, axs
end



function plot_normality_checks(nc_path::String)
    """
    plot_normality_checks(nc_path::String)

    Creates a thesis-ready figure with two plots for multivariate normality diagnostics:
    1. Left: Chi-square PDF vs Empirical Histogram (of squared Mahalanobis distances).
    2. Right: Q-Q Plot (Theoretical vs Sample Quantiles).

    # Arguments
    - `nc_path`: Path to the NetCDF file containing `chi2_x`, `chi2_pdf`, `chi2_hist`, `qq_theor_q`, `qq_D_sorted`.
    """
    # 1. Load Data
    ds = NCDataset(nc_path, "r")

    # Read variables. [:] ensures we get standard Julia vectors
    chi2_x      = ds["chi2_x"][:]
    chi2_pdf    = ds["chi2_pdf"][:]
    chi2_hist   = ds["chi2_hist"][:]
    qq_theor_q  = ds["qq_theor_q"][:]
    qq_D_sorted = ds["qq_D_sorted"][:]

    close(ds)

    # 2. Setup Theme (Consistent with mcmc_plots)
    colors = Makie.to_colormap(:seaborn_dark)
    set_theme!(palette = (; color = colors), fontsize = 24)

    # 3. Initialize Figure
    # Width 1200 to match the side-by-side style of previous plots
    fig = Figure(size = (1100, 400), figure_padding = (5, 20, 5, 5))

    # --- Left Plot: Chi-Square Test ---
    ax_chi = Axis(fig[1, 1],
        # xlabel = L"\mathrm{Squared\ Mahalanobis\ Distance}\ (D^2)",
        # ylabel = L"\mathrm{Density}",
        xgridvisible = true, ygridvisible = true,
        # xticklabelsize = 18, yticklabelsize = 18
    )

    # Plot lines as requested (using explicit colors from your snippet)
    # Slightly thicker lines (2.5) for thesis readability
    lines!(ax_chi, chi2_x, chi2_hist, color = colors[2], linewidth = 4,
           label = L"\mathrm{Empirical\ PDF}")
    lines!(ax_chi, chi2_x, chi2_pdf, color = colors[1], linewidth = 4,
           label = L"\mathrm{Theoretical\ PDF}")

    axislegend(ax_chi, position = :rt, framevisible = false)


    # --- Right Plot: Q-Q Plot ---
    ax_qq = Axis(fig[1, 2],
        # xlabel = L"\mathrm{Theoretical\ Quantiles}\ (\chi^2)",
        # ylabel = L"\mathrm{Sample\ Quantiles}\ (D^2)",
        xgridvisible = true, ygridvisible = true,
        # xticklabelsize = 18, yticklabelsize = 18
    )

    # 1:1 Reference Line (Red dashed)
    ablines!(ax_qq, 0, 1, color = colors[2], linewidth = 4, linestyle = :dash)# label = L"1:1\ \mathrm{Line}")

    # Scatter points (Blue)
    scatter!(ax_qq, qq_theor_q, qq_D_sorted, color = colors[1], markersize = 10,
    )#label = L"\mathrm{Data\ Points}")

    # axislegend(ax_qq, position = :lt, framevisible = false)

    # 4. Final Layout Adjustments
    colgap!(fig.layout, 30)

    return fig, [ax_chi, ax_qq]
end


function plot_model_predictions(nc_path::String, model_name::String)
    """
    plot_model_predictions_with_ci(nc_path::String, model_name::String)

    Generates a plot of model predictions with 95% confidence intervals.
    The model is specified by `model_name` ("OU", "L3", or "blowfly").
    """

    # 1. Load Data
    ds = NCDataset(nc_path, "r")
    chain = ds["chain"][:, :]
    R0_all = ds.attrib["R0_all"]


    # Parse R0_all (vector of vector of matrices)
    R0_all_parsed = eval(Meta.parse(R0_all))
    ground_truth = R0_all_parsed[1]  # First matrix of the vector of vector of matrices

    # 2. Randomly Sample 1000 Parameter Vectors
    n_samples = 1000
    n_chain_samples = size(chain, 1)
    param_indices = rand(1:n_chain_samples, n_samples)
    sampled_params = chain[param_indices, :]

    # 3. Generate Data Based on Model
    model_func = nothing
    realizations = Vector{Matrix{Float64}}(undef, 1000)
    if model_name == "OU"
        model_func = OU_solve
        init = 3.0
        Ndata = 200
        for ii in 1:1000
            params = sampled_params[ii, :]
            realizations[ii] = model_func( params, init, Ndata)
        end
    elseif model_name == "L3"
        # model_func = lorenz_solve

    elseif model_name == "blowfly"
        model_func = blowfly_solve
        init = 180  # Initial condition for synthetic data
        N = 200     # Length of synthetic data
        burn_in = 0    # Time after which the data is considered "stable"
        t = N+1  # Time vector from 0 to N
        for ii in 1:1000
            params = sampled_params[ii, :]
            realizations[ii] = model_func(params, t; n_init = init, burn_in = burn_in)[1]
        end
    else
        error("Invalid model name. Options: 'OU', 'L3', 'blowfly'")
    end
    close(ds)

    if model_name == "OU" || model_name == "blowfly"

        realizations = vcat( realizations... )  # Convert vector of matrices to a single matrix
        # 4. Calculate 95% Confidence Interval
        lower_ci = map(x -> quantile( vec(x), 0.025), eachcol(realizations))
        upper_ci = map(x -> quantile( vec(x), 0.975), eachcol(realizations))
        mean_prediction = mean(realizations, dims=1)

        # 5. Plot Ground Truth and 95% CI
        colors = Makie.to_colormap(:seaborn_dark)
        set_theme!(palette = (; color = colors), fontsize = 21)

        fig = Figure(size = (900, 350), figure_padding = (5, 5, 5, 5))
        ax = Axis(fig[1, 1],
            xlabel = L"\mathrm{Time}",
            ylabel = L"\mathrm{Value}",
            xgridvisible = true, ygridvisible = true,
            xlabelsize = 24, ylabelsize = 24,
        )

        time_points = 1:size(ground_truth, 2)
        lines!(ax, time_points, ground_truth[1, :],
                color = colors[1],
                linewidth = 3,
                label = L"\mathrm{Ground\ Truth}")
        band!(ax, time_points, lower_ci, upper_ci,
                color = (colors[1], 0.3),
                label = L"\mathrm{95\%\ CI}")
        lines!(ax, time_points[:], mean_prediction[:],
                color = colors[3],
                linewidth = 2,
                linestyle = :dash,
                label = L"\mathrm{Mean}")

        axislegend(ax, position = :rt, framevisible = false)

    else
        Ndata = 200+10
        init = 12.577, 19.471, 23.073
        dt = 1.0
        for ii in 1:1000
            params = sampled_params[ii, :]
            realizations[ii] = lorenz_solve( init, params, Ndata; dt = dt )
        end
        # realizations = vcat( realizations... )
        realizations = sum( realizations )
        realizations = realizations ./ 1000.0
        t = ( 0:Ndata ) .* dt
        t = t[11:end]

        fig = Figure(size = (1200, 400), figure_padding = (1, -70, 1, -28), fontsize = 26)
        ax_ts = Axis(fig[1, 1],
            xlabel = L"\mathrm{Time}",
            xgridvisible = true, ygridvisible = true,
            xticklabelsize = 22, yticklabelsize = 22
        )
        # Plot x, y, z time series
        labels_raw = ["x", "y", "z"]
        for i in 1:3
            # Create upright LaTeX label: $\mathrm{x}$, $\mathrm{y}$, etc.
            lbl = latexstring("\\mathrm{" * labels_raw[i] * "}")
            lines!(ax_ts, t, ground_truth[i, :], linewidth = 2, label = lbl)
            lines!(ax_ts, t, R[i, :], linewidth = 2, label = lbl)
        end
        axislegend(ax_ts, position = (0.5, 0.025), framevisible = false)


        # --- 4. Right Panel: Butterfly Attractor (3D) ---
        # Spanning column 3 (1/3rd of width)
        ax_3d = Axis3(fig[1, 2],
            azimuth = -8*pi/19,
            elevation = 0.4, # Adjusted for better view
            xlabel = L"\mathrm{x}",
            ylabel = L"\mathrm{y}",
            zlabel = L"\mathrm{z}",
            xlabeloffset = 30, zlabeloffset = 40, ylabeloffset = 50,
            xticklabelsize = 22, yticklabelsize = 22, zticklabelsize = 22,
            # viewmode = :fit
        )

        # Plot the full trajectory (or a subset if it's too dense)
        lines!(ax_3d, R[1, :], R[2, :], R[3, :], linewidth = 1.0 )
        lines!(ax_3d, ground_truth[1, :], ground_truth[2, :], ground_truth[3, :], linewidth = 1.0, color = :red)

        # --- 5. Final Adjustments ---
        colgap!(fig.layout, 1)
    end

    return fig, ax
end
