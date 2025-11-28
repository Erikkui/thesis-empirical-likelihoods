function fig = mcmcplot( data, NameValueArgs ) %variable, plot_type, labels )
%MCMCPLOT Create plots from NetCDF file 
%   data = array containing names of NetCDF files or mcmc chain matrices
%   variable = desired variable to extract from the file(s)
%   plot_type = type of plot wanted ("chain", "histogram", "pairs")
%   labels = array containing labels for each file in ncfiles

    arguments
        data 
        NameValueArgs.variable
        NameValueArgs.plot_type
        NameValueArgs.labels
        NameValueArgs.burn_in
    end

    ndata = length( data );
    plot_type = NameValueArgs.plot_type;

    if isfield( NameValueArgs, "labels" )
        labels = NameValueArgs.labels;
    else
        labels = string( 1:ndata );
    end

    if isfield( NameValueArgs, "variable" )
        variable = NameValueArgs.variable;
    end

    if isfield( NameValueArgs, "burn_in" )
        burn_in = NameValueArgs.burn_in;
    else
        burn_in = 1;
    end
  
    datas = {};
    for ii = 1:ndata
        obj = data{ii};
        if ischar( obj ) || isstring( obj )
            data_temp = ncread( obj, variable );
            try
                theta = str2num( ncreadatt( data(ii), "/", "params") );
            catch
                theta = data_temp(1, :);
            end
            datas{ii} = data_temp( :, : );
        elseif ismatrix( obj )
            datas{ii} = obj;
        end
    end

    cols = size( datas{1}, 2 );

    ax_cols = 2;
    ax_rows = ceil( cols / ax_cols );
    fig = figure();
    
    switch plot_type
        case "chain"
            tiledlayout( ax_rows, ax_cols );
            inds = 0;
            for dd = 1:ndata
                DD = datas{dd};
                inds = inds(end)+1:inds(end)+size( DD, 1 );

                for ii = 1:cols
                    ax = nexttile;
                    hold on
                    plot( inds, DD( :, ii ), ".", "DisplayName", labels(dd) )
                    ylabel( ['$\theta_{', num2str(ii), '}$'], 'Interpreter', 'latex' )
                    if labels(1) ~= "0"
                        legend(ax, "show")
                    else
                        xlim([0, size(DD, 1)])
                    end
                end
            end

        case "histogram"
            datas = datas{1};
            tiledlayout( ax_rows, ax_cols );
            for ii = 1:cols
                nexttile(ii)
                hold on
                histogram( datas( burn_in:end, ii), 30 );
                xline( theta(ii), "LineWidth", 3, "Color", "red" )
                xlabel( ['$\theta_{', num2str(ii), '}$'], 'Interpreter', 'latex' )
            end

        case  "pairs"
            tiledlayout;
          
            for dd = 1:ndata
                DD = datas{dd};
                for ii = 2:cols
                    ax = nexttile(ii-1);
                    hold on
                    xx = DD( burn_in:end, 1 );
                    yy = DD( burn_in:end, ii );
                    xlabel( '$\theta_{1}$', 'Interpreter', 'latex' )
                    ylabel( ['$\theta_{', num2str(ii), '}$'], 'Interpreter', 'latex' )
                    scatter( xx, yy , "filled", "DisplayName", labels(dd) )
                    plot( theta( 1 ), theta( ii ), 'xr', ...
                             "HandleVisibility", "off", ...
                             "MarkerSize", 8, "LineWidth", 1.5)
                    if labels(1) ~= "0"
                        legend(ax, "show")
                    end
                end
            end

    end

    theme( fig, "light" )

    fh = findobj('-property', 'FontName');
    set( fh, 'FontSize', 24)
    fh = findall(0, 'Type', 'Axes');
    set( fh, 'TitleFontSizeMultiplier', 1.5)
end


