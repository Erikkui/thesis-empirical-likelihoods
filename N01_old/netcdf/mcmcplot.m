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
    end

    ndata = length( data );
    plot_type = NameValueArgs.plot_type;

    if isfield( NameValueArgs, "labels" )
        labels = NameValueArgs.labels;
    else
        labels = string( 1:ndata );
    end

    if iscell(data)
        CASE = "array";
        theta = data{1}(1, :);
        datas = data;
    elseif ischar(data) || isstring(data)
        CASE = "netcdf";
        variable = NameValueArgs.variable;
        datas = {};
        for ii = 1:ndata
            datas{ii} = ncread( data(ii), variable );
        end
    
        try
            theta = str2num( ncreadatt( data(ii), "/", "params") );
        catch
            theta = datas{1}(1, :);
        end
    else
        error('Data must be a cell array of matrices or an array of strings.');
    end
    
    cols = size( datas{1}, 2 );
    fig = figure();
    
    switch plot_type
        case "chain"
            tiledlayout( cols, 1);
            inds = 0;
            for dd = 1:ndata
                DD = datas{dd};
                inds = inds(end)+1:inds(end)+size( DD, 1 );

                for ii = 1:cols
                    ax = nexttile(ii);
                    hold on
                    scatter( inds, DD( :, ii ), "filled", "DisplayName", labels(dd) )
                    ylabel( ['$\theta_{', num2str(ii), '}$'], 'Interpreter', 'latex' )
                    legend(ax, "show")
                end
            end

        case "histogram"
            datas = datas{1};
            tiledlayout( 1, cols );
            for ii = 1:cols
                nexttile(ii)
                hold on
                histogram( datas(:, ii), 20 );
                xline( theta(ii), "LineWidth", 3, "Color", "red" )
                xlabel( theta(ii) )
            end

        case  "pairs"
            tiledlayout;
          
            for dd = 1:ndata
                DD = datas{dd};
                for ii = 2:cols
                    ax = nexttile(ii-1);
                    hold on
                    xx = DD( :, 1 );
                    yy = DD( :, ii );
                    xlabel( '$\theta_{1}$', 'Interpreter', 'latex' )
                    ylabel( ['$\theta_{', num2str(ii), '}$'], 'Interpreter', 'latex' )
                    scatter( xx, yy , "filled", "DisplayName", labels(dd) )
                    plot( theta( 1 ), theta( ii ), 'xr', ...
                             "HandleVisibility", "off", ...
                             "MarkerSize", 8, "LineWidth", 1.5)
                    legend(ax, "show")
                end
            end

    end

    theme( fig, "light" )

    fh = findobj('-property', 'FontName');
    set( fh, 'FontSize', 14)
    fh = findall(0, 'Type', 'Axes');
    set( fh, 'TitleFontSizeMultiplier', 1.5)
end


