clc; close all; clearvars;

% Chain from NetCDF file; specify as needed
chain = ncread( "blowfly_bsl_cov_dt1_n109n1_s10r500_eCDF0_set2.nc", "chain" );
chain_burn = chain( 5000:end, : );  % If need to pick only portion of the chain

% Read data from csv
filepath = "/home/eki/GitHub/likelihood-free-methods/blowfly/materials/blowflies.csv";
data_all = readmatrix( filepath ); 
setnum = 2;     % Which dataset to use
dataset_mask = data_all(:, end) == setnum;
dataset = data_all( dataset_mask, 2 );  

% PLot the dataset
fig = figure();
plot( 0:numel(dataset)-1, dataset, "LineWidth", 3, "Color", "r" )
hold on 

%%
% Pick random parameters from the chain and generate a signal to compare
% with dataset
tsim = numel( dataset );
for ii = 1:5
    theta_ii_ind = randi( size(chain_burn, 1) );
    theta_ii = chain_burn( theta_ii_ind, : );

    [N, ~] = blowfly_solve( theta_ii, tsim );
    tau = floor( theta_ii(5) );
    tt = tau+1:tsim;
    plot( tt-tau-1, N(tt), "-k", "LineWidth", 1)
end



