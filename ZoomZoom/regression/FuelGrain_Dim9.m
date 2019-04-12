close all; clear all; clc

%center hole
x0 = 0;
y0 = 0; 
r0 = 0.125; 

% inner ring
n1 = [1:1:7];
ring1 = 0.39; 
d01 = 45;
x1 =  ring1.*cosd( n1.* d01);
y1 =  ring1.*sind( n1.* d01);
r1 = 0.185/2;

% outer ring
n2 = [1:1:15];
ring2 = 0.725; 
d02 = 22.5;
x2 =  ring2.*cosd( n2.* d02);
y2 =  ring2.*sind( n2.* d02);
r2 = 0.1718/2;

x = [x0, x1, x2]';
y = [y0, y1, y2]';
r = [r0, r1*ones(size(x1)), r2*ones(size(x2))];

filename = 'design1.mat';
save(filename)