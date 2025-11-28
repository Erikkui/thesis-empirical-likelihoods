clc; close all

data = readmatrix( "blowflies.csv" );

setcol = data( :, end );
E1 = data( setcol == 2, 1:2 );
E2 = data( setcol == 3, 1:2 );
E3 = data( setcol == 4, 1:2 );

f = figure;
tiledlayout
nexttile(1)
plot(E1(:, 1), E1(:, 2), 'r');
hold on
nexttile(2)
plot(E2(:, 1), E2(:, 2), 'g');
hold on
nexttile(3)
plot(E3(:, 1), E3(:, 2), 'b');

