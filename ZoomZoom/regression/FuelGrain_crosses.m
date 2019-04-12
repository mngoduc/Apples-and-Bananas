close all; clear all; clc

%center hole
x0 = 0;
y0 = 0; 

r0 = 0.25/2; 

% ring 1 
holes1 = 4;
n1 = [1:1:holes1];
ring1 = 0.25; 
d01 = 360/holes1;
x1 =  ring1.*cosd(45+ n1.* d01);
y1 =  ring1.*sind(45+ n1.* d01);
r1 = 0.1015/2;

% ring 2 
holes2 = 4;
n2 = [1:1:holes2];
ring2 = 0.3; 
d02 = 360/holes2;
x2 =  ring2.*cosd(n2.* d02);
y2 =  ring2.*sind(n2.* d02);
r2 = 0.1562/2;

% ring 3 
holes3 = 4;
n3 = [1:1:holes3];
ring3 = 0.5; 
d03 = 360/holes3;
x3 =  ring3.*cosd(45+n3.* d03);
y3 =  ring3.*sind(45+n3.* d03);
r3 = 0.1890/2;

% ring 4 
holes4 = 4;
n4 = [1:1:holes4];
ring4 = ring2 + r2 + 0.15; 
d04 = 360/holes4;
x4 =  ring4.*cosd(n4.* d04);
y4 =  ring4.*sind(n4.* d04);
r4 = 0.1470/2;

% ring 5 
holes5 = 4;
n5 = [1:1:holes5];
ring5 = ring4 + r4 + 0.15; 
d05 = 360/holes5;
x5 =  ring5.*cosd(n5.* d05);
y5 =  ring5.*sind(n5.* d05);
r5 = 0.120/2;

% ring 6 
holes6 = 4;
n6 = [1:1:holes6];
ring6 = 0.75; 
d06 = 360/holes6;
x6 =  ring6.*cosd(-25+ n6.* d06);
y6 =  ring6.*sind(-25+ n6.* d06);
r6 = 0.1015/2;

% ring 7 
holes7 = 4;
n7 = [1:1:holes7];
ring7 = 0.75; 
d07 = 360/holes7;
x7 =  ring7.*cosd(25+ n7.* d07);
y7 =  ring7.*sind(25+ n7.* d07);
r7 = 0.1015/2;

x = [x0, x2, x3, x4, x5, x6, x7]';
y = [y0, y2, y3, y4, y5, y6, y7]';
r = [r0,...
    r2*ones(size(x2)), ...
    r3*ones(size(x3)),...
    r4*ones(size(x4)),...
    r5*ones(size(x5)),...
    r6*ones(size(x6)),...
    r7*ones(size(x7))];

filename = 'XOXO.mat';
save(filename)