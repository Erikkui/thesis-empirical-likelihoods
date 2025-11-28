function create_figures(chain, theta, fig_path; show=true, savefig=false, case_name="test")
    AX_COLS = 3

    n_chain = size( chain, 1 )  # Length of the chain
    npar = size(chain, 2)  # Number of parameters

    chain_burned = chain[ Int(round(n_chain/3)):end, : ]  # Burn-in third the chain
    ax_rows = max( 1, Int( ceil( (npar-1)/AX_COLS ) ) )
    with_theme(theme_latexfonts()) do
        ### Posteriors
        fig2 = Figure( size = ( 1280, 400*ax_rows ) )

        fig2_axes = Vector{Axis}( undef, npar-1 )
        for ii in 1:ax_rows
            for jj in 1:AX_COLS
                ax_ind = (ii-1)*AX_COLS + jj
                theta_ind = string( ax_ind+1 )
                xlabel = L"\theta_1"
                ylabel = L"\theta_%$(theta_ind)"
                ax = Axis( fig2[ii, jj], xlabel = xlabel, ylabel = ylabel,
                    xlabelsize = 30, ylabelsize = 30, xticklabelsize = 25, yticklabelsize = 25 )

                scatter!( ax, chain_burned[:, 1], chain_burned[:, ax_ind+1] )
                scatter!( ax, theta[1], theta[ax_ind+1], markersize = 20, marker = :xcross, color = :red )

                fig2_axes[ ax_ind ] = ax

                if ax_ind == npar-1
                    break
                end
            end
        end
        fig_title = "Parameter posteriors for " * CASE_NAME * " model"
        fig2[0, :] = Label( fig2, fig_title, fontsize = 40 )

        fig_name = string( "run_", case_name, "_posteriors.pdf" )
        save_path = joinpath( fig_path, fig_name )
        if savefig
            save( save_path, fig2 )
        end
        if show
            display( fig2 )
        end


        ### Marginals
        fig3 = Figure( size = (1280, 400*ax_rows) )
        fig3_axes = Vector{Axis}( undef, npar )
        ax_cols = min( npar, AX_COLS )
        for ii in 1:ax_rows
            for jj in 1:ax_cols
                ax_ind = (ii-1)*ax_cols + jj
                theta_ind = string( ax_ind )
                xlabel = L"\theta_%$(theta_ind)"
                ax = Axis( fig3[ii, jj], xlabel = xlabel,
                    xlabelsize = 30, ylabelsize = 30, xticklabelsize = 20, yticklabelsize = 20 )

                hist!( ax, chain_burned[:, ax_ind], bins = 15, color = :blue, strokewidth = 1, strokecolor = :black )
                vlines!( ax, [ theta[ax_ind] ], color = :red, linewidth = 3 )
                fig3_axes[ ax_ind ] = ax
            end
        end

        fig_name = string( "run_", case_name, "_marginals.pdf" )
        save_path = joinpath( fig_path, fig_name )
        if savefig
            save( save_path, fig3 )
        end
        if show
            display( fig3 )
        end


        ### Chains
        fig4 = Figure(size = (1000, npar*200))
        fig4_axes = Vector{Axis}(undef, npar)
        for ii in 1:npar
            ii_string = string(ii)
            ax = Axis( fig4[ii, :], ylabel = L"\theta_%$(ii_string)",
                xlabelsize = 20, ylabelsize = 20, xticklabelsize = 15, yticklabelsize = 15,
                width = Relative(3.0), height = Relative(1.0) )

            scatter!(ax, 1:size(chain, 1), chain[:, ii], color = :blue )
            fig4_axes[ii] = ax
        end
        fig4[0, :] = Label(fig4, "MCMC chains", fontsize = 30)

        fig_name = string( "run_", case_name, "_chains.pdf" )
        save_path = joinpath( fig_path, fig_name )
        if savefig
            save( save_path, fig4 )
        end
        if show
            display( fig4 )
        end
    end

    return nothing
end
