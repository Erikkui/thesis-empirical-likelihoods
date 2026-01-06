using CairoMakie
using NCDatasets
using LaTeXStrings
using MathTeXEngine
using Statistics


function plot_mcmc_results(nc_path::String, plot_type::Symbol, param_names::Vector{String}, true_params; burn_in::Int=0)
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
    close(ds)

    # 2. Apply Burn-in
    if burn_in >= size(raw_chain, 1)
        println("Burn-in length ($burn_in) exceeds chain length.")
        burn_in = 0
    end
    chain = raw_chain[burn_in+1:end, :]
    n_total_samples, n_params = size(raw_chain)

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
        full_vec = raw_chain[:, i]
        burned_vec = chain[:, i]
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
                yticks = Makie.WilkinsonTicks(4), # Same tick logic as chains
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
            lines!(ax_chain, 1:n_total_samples, full_vec, linewidth=1.5)

            # Burn-in Line
            if burn_in > 0
                vlines!(ax_chain, [burn_in], color = (:grey, 0.8), linestyle = :dash, linewidth = 2.0)
            end

            # Rotated Density using burned chain (direction=:x makes it horizontal)
            hist!(ax_hist, burned_vec, normalization = :pdf, direction = :x, strokewidth = 1, color = colors[3], strokecolor = :black, bins = 20)

            # True Parameter Lines (Horizontal)
            if !(contains( nc_path, "set" ))
                hlines!(ax_chain, [true_params[i]], color = (colors[2], 0.85), linewidth = 3.75)
                hlines!(ax_hist, [true_params[i]], color = (colors[2], 0.85), linewidth = 3.75)
            end

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
    fig = Figure(size = (1100, 400), figure_padding = (5, 5, 5, 10))

    # --- Left Plot: Chi-Square Test ---
    ax_chi = Axis(fig[1, 1],
        # xlabel = L"\mathrm{Squared\ Mahalanobis\ Distance}\ (D^2)",
        # ylabel = L"\mathrm{Density}",
        xgridvisible = true, ygridvisible = true,
        # xticklabelsize = 18, yticklabelsize = 18
    )

    # Plot lines as requested (using explicit colors from your snippet)
    # Slightly thicker lines (2.5) for thesis readability
    lines!(ax_chi, chi2_x, chi2_hist, color = colors[2], linewidth = 4)
        #    label = L"\mathrm{Empirical\ PDF}")
    lines!(ax_chi, chi2_x, chi2_pdf, color = colors[1], linewidth = 4)
        #    label = L"\mathrm{Theoretical\ PDF}")

    # axislegend(ax_chi, position = :rt, framevisible = false)


    # --- Right Plot: Q-Q Plot ---
    ax_qq = Axis(fig[1, 2],
        # xlabel = L"\mathrm{Theoretical\ Quantiles}\ (\chi^2)",
        # ylabel = L"\mathrm{Sample\ Quantiles}\ (D^2)",
        xgridvisible = true, ygridvisible = true,
        # xticklabelsize = 18, yticklabelsize = 18
    )

    # 1:1 Reference Line (Red dashed)
    ablines!(ax_qq, 0, 1, color = colors[1], linewidth = 4, linestyle = :dash)# label = L"1:1\ \mathrm{Line}")

    # Scatter points (Blue)
    scatter!(ax_qq, qq_theor_q, qq_D_sorted, color = colors[2], markersize = 10,
    )#label = L"\mathrm{Data\ Points}")

    # axislegend(ax_qq, position = :lt, framevisible = false)

    # 4. Final Layout Adjustments
    colgap!(fig.layout, 30)

    return fig, [ax_chi, ax_qq]
end


function plot_model_predictions(nc_path::String, model_name::String; burn_in::Int=0, n_samples = 100)
    """
    plot_model_predictions_with_ci(nc_path::String, model_name::String)

    Generates a plot of model predictions with 95% confidence intervals.
    The model is specified by `model_name` ("OU", "L3", or "blowfly").
    """

    # 1. Load Data
    ds = NCDataset(nc_path, "r")
    chain = ds["chain"][:, :]
    R0_all = ds.attrib["R0_all"]

    if burn_in >= size(chain, 1)
        println("Burn-in length ($burn_in) exceeds chain length.")
        burn_in = 0
    end
    chain = chain[burn_in+1:end, :]

    # Parse R0_all (vector of vector of matrices)
    R0_all_parsed = eval(Meta.parse(R0_all))
    ground_truth = R0_all_parsed[1]  # First matrix of the vector of vector of matrices

    # 2. Randomly Sample 1000 Parameter Vectors + posterior mean
    n_chain_samples = size(chain, 1)
    param_indices = rand(1:n_chain_samples, n_samples)
    sampled_params = chain[param_indices, :]
    posterior_mean = mean(chain, dims=1)
    sampled_params = vcat( sampled_params, posterior_mean )  # Add posterior mean as
    n_samples += 1  # Update sample count

    # 3. Generate Data Based on Model
    model_func = nothing
    realizations = Vector{Matrix{Float64}}(undef, n_samples)
    if model_name == "OU"
        init = 3.0
        Ndata = 200
        for ii in 1:n_samples
            params = sampled_params[ii, :]
            realizations[ii] = OU_solve( params, init, Ndata)
        end
    elseif model_name == "L3"
        # model_func = lorenz_solve

    elseif model_name == "blowfly"
        init = parse(Float64, ds.attrib["synth_init"])  # Initial condition for synthetic data
        t = parse(Int, ds.attrib["synth_N"])     # Length of synthetic data
        burn_in = 0    # Time after which the data is considered "stable"
        # t = N+1  # Time vector from 0 to N
        for ii in 1:n_samples
            params = sampled_params[ii, :]
            realizations[ii] = blowfly_solve(params, t; N_init = init, burn_in = burn_in)[1]
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
        set_theme!(palette = (; color = colors), fontsize = 20)

        fig = Figure(size = (900, 350), figure_padding = (5, 5, 5, 5))
        ax = Axis(fig[1, 1],
            xlabel = L"\mathrm{Time}",
            ylabel = L"\mathrm{Value}",
            xgridvisible = true, ygridvisible = true,
            xlabelsize = 23, ylabelsize = 23,
        )

        time_points = 1:size(ground_truth, 2)
        lines!(ax, time_points, ground_truth[1, :],
                color = colors[1],
                linewidth = 3.5,
                label = L"\mathrm{Ground\ Truth}")
        if model_name == "blowfly"
            lines!(ax, time_points[:], realizations[ end, : ],
            color = colors[2],
            linewidth = 2.5,
            linestyle = :solid,
            label = L"\mathrm{Posterior\ mean}")
        end
        # lines!(ax, time_points[:], mean_prediction[:],
        #         color = colors[3],
        #         linewidth = 2,
        #         linestyle = :dash,
        #         label = L"\mathrm{Mean}")
        band!(ax, time_points, lower_ci, upper_ci,
                color = (colors[3], 0.2),
                label = L"\mathrm{95\% \ CI}")

        Legend(fig[2, 1], ax, orientation = :horizontal, framevisible = false, labelsize=23)
        # colgap!(fig.layout, 0)
    else
        # Ndata = 200+10
        # init = [12.577, 19.471, 23.073]
        # dt = 1.0
        # for ii in 1:n_samples
        #     params = sampled_params[ii, :]
        #     realizations[ii] = lorenz_solve( init, params, Ndata; dt = dt )
        # end
        # x = [ r[1, :] for r in realizations ]
        # x = vcat( x... )
        # y = [ r[2, :] for r in realizations ]
        # y = vcat( y... )
        # z = [ r[3, :] for r in realizations ]
        # z = vcat( z... )
        # R = [ mean(x, dims = 1)'; mean(y, dims = 1)'; mean(z, dims = 1)' ]
        # realizations = vcat( realizations... )
        # realizations = sum( realizations )
        # realizations = realizations ./ n_samples
        # t = ( 0:Ndata ) .* dt
        # t = t[11:end]

        # fig = Figure(size = (1200, 400), figure_padding = (1, -70, 1, -28), fontsize = 26)
        # ax_ts = Axis(fig[1, 1],
        #     xlabel = L"\mathrm{Time}",
        #     xgridvisible = true, ygridvisible = true,
        #     xticklabelsize = 22, yticklabelsize = 22
        # )
        # Plot x, y, z time series
        # labels_raw = ["x", "y", "z"]
        # for i in 1:3
        #     Create upright LaTeX label: $\mathrm{x}$, $\mathrm{y}$, etc.
        #     lbl = latexstring("\\mathrm{" * labels_raw[i] * "}")
        #     lines!(ax_ts, t, ground_truth[i, :], linewidth = 2, label = lbl)
        #     lines!(ax_ts, t, R[i, :], linewidth = 2, label = lbl)
        # end
        # axislegend(ax_ts, position = (0.5, 0.025), framevisible = false)


        # --- 4. Right Panel: Butterfly Attractor (3D) ---
        # Spanning column 3 (1/3rd of width)
        # ax_3d = Axis3(fig[1, 2],
        #     azimuth = -8*pi/19,
        #     elevation = 0.4, # Adjusted for better view
        #     xlabel = L"\mathrm{x}",
        #     ylabel = L"\mathrm{y}",
        #     zlabel = L"\mathrm{z}",
        #     xlabeloffset = 30, zlabeloffset = 40, ylabeloffset = 50,
        #     xticklabelsize = 22, yticklabelsize = 22, zticklabelsize = 22,
        #     viewmode = :fit
        # )

        # Plot the full trajectory (or a subset if it's too dense)
        # lines!(ax_3d, R[1, :], R[2, :], R[3, :], linewidth = 1.0 )
        # lines!(ax_3d, ground_truth[1, :], ground_truth[2, :], ground_truth[3, :], linewidth = 1.0, color = :red)

        # --- 5. Final Adjustments ---
        # colgap!(fig.layout, 1)
    end

    return fig
end






# --- Helper: HPD Calculation ---
function calculate_hpd(samples::Vector{T}, params_true; alpha=0.05) where T <: Real
    n = length(samples)
    m = round(Int, n * (1 - alpha))
    y = sort(samples)
    ranges = y[m+1:n] - y[1:n-m]
    min_idx = argmin(ranges)

    posterior_mean = mean(samples)
    z_score = abs( (params_true - mean(samples)) / std(samples) )

    return (lower = y[min_idx], upper = y[min_idx + m]), posterior_mean, z_score
end


function plot_forest_multi(file_paths::Vector{String}, param_names::Vector{String}, params_true;
                           burn_in::Int=1000,
                           labels::Vector{String}=String[],
                           var_name::String="chain")

    """
        plot_forest_multi(file_paths::Vector{String}, param_names::Vector{String};
                        burnin::Int=1000,
                        labels::Vector{String}=String[],
                        var_name::String="chain")

    Creates a thesis-ready forest plot comparing multiple parameters side-by-side.
    Matches the style of `plot_normality_checks`.

    # Arguments
    - `file_paths`: Vector of paths to NetCDF files.
    - `param_names`: Vector of LaTeX-formatted strings for x-axis labels (e.g., [L"\\theta_1", L"\\theta_2"]).
    - `labels`: Vector of labels for the y-axis (experiment names).
    """

    # 1. Setup Theme (Matching your reference)
    colors = Makie.to_colormap(:seaborn_dark)
    # Use the same blue (1) and red/orange (2) as your reference
    color_main = colors[1]
    color_ref  = colors[2]

    set_theme!(palette = (; color = colors), fontsize = 24)

    # 2. Data Loading & Processing
    num_params = length(param_names)

    # Structure: plot_data[param_index] = Vector of (mean, lower, upper, truth, label)
    plot_data = [ [] for _ in 1:num_params ]

    # We capture the "global truth" for each parameter from the first file
    global_truths = fill(NaN, num_params)

    for (file_idx, path) in enumerate(file_paths)
        println("------------------- Processing file: ", split(path, "/")[end], " -------------------")
        if !isfile(path)
            @warn "File not found: $path"
            continue
        end

        hpds = []
        post_means = []
        z_scores = []

        NCDataset(path, "r") do ds
            if !haskey(ds, var_name)
                error("Variable '$var_name' not found in $path")
            end

            raw_chain = ds[var_name][:, :]
            # Handle 1D vs 2D chains
            chain_matrix = ndims(raw_chain) == 1 ? reshape(raw_chain, :, 1) : raw_chain

            if size(chain_matrix, 2) < num_params
                error("File $path has $(size(chain_matrix, 2)) parameters, but $(num_params) names provided.")
            end

            # Loop over each parameter requested
            for p in 1:num_params
                chain_vec = chain_matrix[:, p]

                # Extract Truth (Step 1)
                # Extract Truth from attribute if available, otherwise use first chain value
                current_true = params_true[p]

                if file_idx == 1
                    global_truths[p] = current_true
                end

                # Burn-in
                if length(chain_vec) <= burn_in
                    println("Burn-in length ($burn_in) exceeds chain length.")
                    burn_in = 0
                end
                posterior = chain_vec[burn_in+1:end]

                # Calculate Stats (HPD)
                p_mean = mean(posterior)
                hpd_int, posterior_mean, z_score = calculate_hpd(posterior, params_true[p])
                push!(hpds, hpd_int)
                push!(post_means, posterior_mean)
                push!(z_scores, z_score)

                lbl = isempty(labels) ? basename(path) : labels[file_idx]

                push!(plot_data[p], (
                    label = lbl,
                    mean  = p_mean,
                    lower = hpd_int.lower,
                    upper = hpd_int.upper,
                    zscore = z_score,
                ))
            end

            println("\n--- Posterior Means ---")
            for (i, pm) in enumerate(post_means)
                println("\$", round(pm, sigdigits=3), "\$")
            end

            println("\n--- Z-Scores ---")
            for (i, z) in enumerate(z_scores)
                println("\$", round(z, sigdigits=3), "\$")
            end

            println("\n--- HPD Intervals ---")
            for (i, hpd) in enumerate(hpds)
                println("\$", round(hpd.lower, sigdigits=3), "\$")
                println("\$", round(hpd.upper, sigdigits=3), "\$")
            end

        end
    end

    # 3. Figure Initialization
    # Adjust width based on number of parameters (e.g., 500 per param)
    max_width = 1100
    max_cols = 3
    n_cols = min(num_params, max_cols)
    n_rows = ceil(Int, num_params / max_cols)
    fig_width = max_width
    fig_height = 400 * n_rows
    fig = Figure(size = (fig_width, fig_height), figure_padding = (5, 20, 5, 5))

    axes = Axis[]

    # 4. Create Subplots
    for p in 1:num_params
        data = plot_data[p]
        n_items = length(data)
        y_pos = 1:n_items

        # Extract vectors
        means = [d.mean for d in data]
        lowers = [d.lower for d in data]
        uppers = [d.upper for d in data]

        # Axis Setup
        # Only show y-tick labels (Experiment Names) on the FIRST plot (Leftmost)
        show_ylabels = (mod(p - 1, max_cols) == 0)
        ax_row = div(p - 1, max_cols) + 1
        ax_col = mod(p - 1, max_cols) + 1

        ax = Axis(fig[ax_row, ax_col],
            xlabel = latexstring("\\mathrm{" * param_names[p] * "}"),
            xticks = Makie.WilkinsonTicks(5),
            yticks = (y_pos, [d.label for d in data]),
            yticksvisible = false,
            ylabel = show_ylabels ? latexstring("\\mathrm{Experiment}") : "",
            yticklabelsvisible = show_ylabels,
            xgridvisible = true,
            ygridvisible = false, # Forest plots usually don't have horizontal grids
            xlabelsize = 28, ylabelsize = 28,
            xminorgridvisible = true
        )
        push!(axes, ax)

        # Reverse Y-axis so first item is at top
        ax.yreversed = true

        # A. Reference Line (True Parameter) - Matching your dashed style (color[2])
        vlines!(ax, [global_truths[p]], color = color_ref, linewidth = 4, linestyle = :dash)

        # B. Interval Bars (The Forest)
        rangebars!(ax, y_pos, lowers, uppers, direction = :x, color = (colors[8], 0.7), linewidth = 4)

        # C. Means (Dots) - Color Coded for Bias
        # Logic: Blue (Standard). Red if truth farther than 2 SDs from mean
        dot_colors = [ (abs(d.zscore) > 2) ? colors[4] : color_main for d in data ]

        scatter!(ax, means, y_pos, color = dot_colors, markersize = 15)
    end

    # 5. Final Layout
    colgap!(fig.layout, 20)

    return fig
end
