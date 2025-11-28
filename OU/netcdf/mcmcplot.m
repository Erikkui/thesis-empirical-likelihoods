function fig = mcmcplot( ncfiles, variable, plot_type, labels )
%MCMCPLOT Create plots from NetCDF file 
%   ncfiles = array containing names of NetCDF files
%   variable = desired variable to extract from the file(s)
%   plot_type = type of plot wanted ("chain", "histogram", "pairs")
%   labels = array containing labels for each file in ncfiles

    nfiles = length( ncfiles );
    data = {};
    for ii = 1:nfiles
        data{ii} = ncread( ncfiles(ii), variable );
    end
    cols = size( data{1}, 2 );

    try
        theta = str2num( ncreadatt( ncfiles(ii), "/", "params") );
    catch
        theta = [ 0.16 6.5 400 0.1 14 0.1 ];
    end

    fig = figure();
    
    switch plot_type
        case "chain"
            tiledlayout( cols, 1);
            inds = 0;
            for dd = 1:nfiles
                DD = data{dd};
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
            data = data{1};
            tiledlayout( 1, cols );
            for ii = 1:cols
                nexttile(ii)
                hold on
                histogram( data(:, ii), 15 );
                xline( theta(ii), "LineWidth", 3, "Color", "red" )
                xlabel( theta(ii) )
            end

        case  "pairs"
            tiledlayout;
          
            for dd = 1:nfiles
                DD = data{dd};
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

    theme(fig,"light")

    fh = findobj('-property', 'FontName');
    set( fh, 'FontSize', 14)
    fh = findall(0, 'Type', 'Axes');
    set( fh, 'TitleFontSizeMultiplier', 1.5)
end