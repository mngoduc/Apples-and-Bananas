function [r m_fuel] = solveradius(mdot_O2,r_init,P,time)
%Regression solver for HDPE
%Single circular port, multiply by N to get total mass burned
%Assume constant mdot_O2
%use equation r=a(G)^n*P^b where G=(mdot_O2/Aport)
% Input(s) - mdot_O2 (kg/s), r_init(inches), P(MPa), time(s)
% Output(s) - r(inches) , fuel_burned(g)

%regression rate constants (HDPE) - TU Delft
a=6.3E-3;
n=0.36;
b=0.22;
%n=0.6;
%b=0;

%convert inches to cm
r_init = r_init * 2.54;
p_HDPE = 950; %kg/m^3
FG_length = 7.350 * (2.54/100); %m

%init arrays
r = zeros(size(time));
dr = zeros(size(time));
Aport = zeros(size(time));

r(1) = r_init;   %cm
Aport_init = pi*(r_init/100)^2;  %m^2
Aport(1) = Aport_init;
dT = time(2)-time(1);

for i=1:length(time)
    G = (mdot_O2(i)/Aport(i)); % kg/(m^2s)
    dr(i) = a*(G^n)*(P(i)^b)*dT;  %cm
    r(i+1) = r(i)+dr(i); %cm
    Aport(i+1) = pi*(r(i+1)/100)^2;  %m^2
end

Aport(end) = [];
r(end)=[];

A_burned = Aport - Aport_init;  %m^3
m_fuel = A_burned * FG_length * p_HDPE;   
m_fuel = m_fuel*1000; %g

r = r./2.54;




