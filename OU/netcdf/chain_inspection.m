clc
% Specify any number of ncfiles, put them in one vector and give legends
% for each file. mcmcplot documnetation can be found from the function file

ncfile1 = "blowfly_v2_gsl_cov_dt1_n151n5_s1r30_eCDF-2-10_synth" + ".nc";
ncfile2 = "blowfly_v2_bsl_cov_dt1_n151n5_s10r20_eCDF-2-10_synth" + ".nc";
ncfile3 = "blowfly_v2_bsl_cov_dt1_n151n1_s10r40_eCDF-2-10_synth" + ".nc";
% ncfile4 = "blowfly_bsl_cov_dt1_n359n1_s5r20_eCDF0_set4" + ".nc";

ncdisp( ncfile1 )
ncdisp( ncfile2 )
ncdisp( ncfile3 )


%%

ncfiles = [ ncfile1, ncfile2, ncfile3];%, ncfile4 ];
legends = ["set 1", "set 2", "set3"]; 

mcmcplot( ncfiles, "chain", "chain", legends );

mcmcplot( ncfiles, "chain", "pairs", legends );
