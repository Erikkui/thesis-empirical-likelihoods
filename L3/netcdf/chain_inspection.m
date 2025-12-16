clc
% Specify any number of ncfiles, put them in one vector and give legends
% for each file. mcmcplot documnetation can be found from the function file

ncfile1 = "OU_bsl_cov_n201n1_s200_r1_dt1_eCDF-1_synth" + ".nc";
% ncfile2 = "blowfly_gsl_cov_n251n50_s100_r1_dt1_eCDF-20_synth" + ".nc";
% ncfile3 = "blowfly_v2_bsl_cov_dt1_n151n1_s10r40_eCDF-2-10_synth" + ".nc";
% ncfile4 = "blowfly_bsl_cov_dt1_n359n1_s5r20_eCDF0_set4" + ".nc";

ncdisp( ncfile1 )
% ncdisp( ncfile2 )
% ncdisp( ncfile3 )


%%

% load( "blowfly_chain2_100nepo350nobs_gsl.mat", "chain2")

datas = {ncfile1};
legends = ["0"];%, "julia"];%, "set 2", "set3"]; 

burn = 30000;

mcmcplot( datas, variable="chain", plot_type="chain", labels=legends, burn_in=burn );

mcmcplot( datas, variable="chain", plot_type="pairs", labels=legends, burn_in=burn );

mcmcplot( datas, variable="chain", plot_type="histogram", labels=legends, burn_in=burn );


%%
setnum = 4;

data0 = str2num(ncreadatt( ncfile1, "/", "R0_all" ));
data = data0(1:201);
N_init = str2num(ncreadatt( ncfile1, "/", "R0_full" ));
N_init = N_init(1);
%%
chain = ncread( ncfile1, "chain" );
% chain = chain( burn:end, : );
%%
tt = 0:length(data)-1;
fig = figure;
hold on

posterior_mean = mean(chain, 1);

t = length(data);

styles = ["-", "--", "-."];
nlines = 500;
linedata = [];
for ii = 1:nlines
    [N, N_burned] = blowfly_solve( posterior_mean, t, burn_in = 0, N_init=N_init );
    linedata(end+1, :) = N_burned;
end

meanline = mean( linedata, 1 );

plot( tt, linedata', Color=[0.678, 0.847, 0.902], LineWidth=2, HandleVisibility="off")
plot( tt, meanline, Color = [0,0,1], Linewidth=3, DisplayName="Mean of simulated trajectories")

legend show

% xFill = [ tt(:); flipud( tt(:)) ];
% yFill = [ N_burned_up(:); flipud( N_burned_low(:))];
% fill(xFill, yFill, 'b', "FaceAlpha", 0.2);
xlim([0, length(data)-1])

plot( tt, data, 'r-', LineWidth=4, DisplayName="Original data" )

theme( fig, "light" )
fh = findobj('-property', 'FontName');
set( fh, 'FontSize', 24)
fh = findall(0, 'Type', 'Axes');
set( fh, 'TitleFontSizeMultiplier', 1.5)
