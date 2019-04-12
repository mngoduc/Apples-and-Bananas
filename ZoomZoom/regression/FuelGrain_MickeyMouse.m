close all; clear all; clc

%center hole
x0 = 0;
y0 = 0; 
r0 = 0.3/2; 

% ring 1 
n1 = [1:1:5];
ring1 = 0.30; 
d01 = 360/5;
x1 =  ring1.*cosd(9+18+ n1.* d01);
y1 =  ring1.*sind(9+18+ n1.* d01);
r1 = 0.152/2;

% ring 2 
n2 = [1:1:12];
ring2 = 0.45; 
d02 = 360/12;
x2 =  ring2.*cosd( 18 + n2.* d02);
y2 =  ring2.*sind( 18 + n2.* d02);
r2 = 0.136/2;

% ring 3
n3 = [1:1:10];
ring3 = 0.725; 
d03 = 360/10;
x3 =  ring3.*cosd( -18 + n3.* d03);
y3 =  ring3.*sind( -18 + n3.* d03);
r3 = 0.2/2;

x = [x0, x1, x2, x3]';
y = [y0, y1, y2, y3]';
r = [r0, r1*ones(size(x1)), r2*ones(size(x2)), r3*ones(size(x3))];

filename = 'MickeyMouse.mat';
save(filename)