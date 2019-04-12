%Get better estimation of "n" coefficient using trial 1 data
%Assume a=6.3E-3 and b=0.22 [from TU Delft datasheet]

clear all; close all; clc;

a=6.3E-3;
b=0.22;
i=1;
m_burned_lab = 192.9-146.6;  %g
fuelGrain = 'FuelGrain_Dim6';
labData = 'hotplay2.mat';

for n=0.35:0.01:0.45
    n_coeff(i) = n;
    m_fuel(i) = get_fuel_burned(a,n,b, fuelGrain ,labData);
    i=i+1;
end

m_lab_line = m_burned_lab .* ones(size(n_coeff));

p = polyfit(n_coeff,m_fuel,1);
m_fuel_fit = polyval(p,n_coeff);
m_fuel_fit(m_fuel_fit<0)=0;

plot(n_coeff,m_fuel,'k-.',n_coeff,m_fuel_fit,'k-',n_coeff,m_lab_line,'r','LineWidth',1.1);

%results -> a= 6.3E-3, n= ~0.4, b= 0.22