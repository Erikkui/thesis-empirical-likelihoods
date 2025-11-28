using CairoMakie
using NCDatasets

file1 = "/home/eki/GitHub/likelihood-free-methods/blowfly/netcdf/blowfly_bsl_cov_dt1_n109n1_s10r500_eCDF0_set2.nc"
file2 = "blowfly/netcdf/blowfly_bsl_cov_dt1_n109n1_s10r500_eCDF0_set3.nc"
file3 = "blowfly/netcdf/blowfly_bsl_cov_dt1_n359n1_s5r20_eCDF0_set4.nc"

files = [file1, file2, file3]

npar = size( NCDataset(file1)["chain"], 2 )
fig = Figure( size = (2000, 500 ) )
colors = Makie.wong_colors()[1:length(files)]
axs = [ Axis(fig[1, ii]) for ii in 1:npar-1 ]

labels = ["set 2", "set 3", "set 4"]

for (jj, file) in enumerate( files )
    ds = NCDataset(file)
    local chain = ds["chain"][:, :]
    for ii in 1:npar-1
        ax = axs[ii]
        ax.xlabel = "theta1"
        ax.ylabel = "theta" * string(ii + 1)
        ax.xticklabelsize = 20
        ax.yticklabelsize = 20
        ax.xlabelsize = 20
        ax.ylabelsize = 20
        scatter!(ax, chain[:, 1], chain[:, ii+1], color = colors[jj], label = labels[jj])
        if jj == length(files)
            axislegend(ax)
        end
    end
end

display(fig)

##

fig2 = Figure( size = (1000, 1200 ) )
axs2 = [ Axis(fig2[ii, 1]) for ii in 1:npar ]
for ii in 1:npar
    local chain = NCDataset( files[2] )["chain"][:, :]
    ax = axs2[ii]
    ax.ylabel = "theta_"*string(ii)
    ax.xticklabelsize = 20
    ax.yticklabelsize = 20
    ax.xlabelsize = 20
    ax.ylabelsize = 20
    scatter!( ax, 1:size(chain, 1), chain[:, ii], color = colors[2] )
end

display(fig2)

##
chainn = NCDataset( files[2] )["chain"][:, :]
