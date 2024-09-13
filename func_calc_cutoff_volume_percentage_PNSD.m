function [percent_cutoff_PM1, percent_cutoff_PM25, percent_cutoff_PM10]=func_calc_cutoff_volume_percentage_PNSD(cutoff_D, accum_DG, accum_theta, accum_V, coarse_DG, coarse_theta, coarse_V)
% Copyright (C) Ying Chen, cycy1320@gmail.com

% INPUTS:
% cutoff_D= dry PM diameter [nm] that larger particles will be blocked outside by sampler; e.g., cutoff_D for PM1 = 1000/growth_factor, because dry_PM_Diameter * growth_factor = 1000 nm (cutoff threshold)
% accum_DG = Geometric mean diameter for accumulation mode 
% accum_theta = standard deviation for accumulation mode 
% accum_V = volume concentration for accumulation mode [um3/cm3];  this is no important if calculate for percentage of cutoff, can set as '1'
% coarse_DG, coarse_theta, coarse_V -- similar as accumulation mode, but parameters for coarse model particles
% Note: aiken mode is note considered here, because Aiken mode is not expected to grow to large enough to be blocked out.


% OUTPUTS:
% Hygroscopicity-induced extra cutoff [%]:  percent_cutoff_PM1, percent_cutoff_PM25, percent_cutoff_PM10


dlogDp=0.01;

%% Typical Average Background:  Whitby, 1997, AE
% accum_DG=320;  % unit: nm
% accum_theta=2.0;
% accum_V=4.45; % unit: um3/cm3
% coarse_DG=6040;
% coarse_theta=2.16;
% coarse_V=25.9; % unit: um3/cm3



x1=log10(10):dlogDp:log10(1000);
X=10.^(x1+dlogDp/2);

% accum
for i=1:numel(x1)
    y_accum(i)=X(i)^3 * accum_V*(1/sqrt(2*pi)/log(accum_theta))*exp(-1/2/(log(accum_theta)^2)*(log(X(i))-log(accum_DG))^2);
end
% % coarse
% for i=1:numel(x1)
%     y_coarse(i)=X(i)^3 * coarse_V*(1/sqrt(2*pi)/log(coarse_theta))*exp(-1/2/(log(coarse_theta)^2)*(log(X(i))-log(coarse_DG))^2);
% end
% Y=y_coarse+y_accum;
Y=y_accum;


dV=Y*dlogDp;

% volume in PM1
temp=abs(X-(1000+1e5));
M_PM1=find(temp==min(temp));
dV_PM1=sum(dV(1:M_PM1));

% volume in PM2.5
temp=abs(X-(2500+1e5));
M_PM25=find(temp==min(temp));
dV_PM25=sum(dV(1:M_PM25));

% volume in PM10
temp=abs(X-(10000+1e5));
M_PM10=find(temp==min(temp));
dV_PM10=sum(dV(1:M_PM10));




temp=abs(X-cutoff_D);
M_cutoff=find(temp==min(temp));
%% cutoff volume in PM1
    dV_cutoff_PM1=sum(dV(M_cutoff:M_PM1));
    percent_cutoff_PM1=dV_cutoff_PM1/dV_PM1*100; %unit: %
%% cutoff volume in PM2.5
    dV_cutoff_PM25=sum(dV(M_cutoff:M_PM25));
    percent_cutoff_PM25=dV_cutoff_PM25/dV_PM25*100; %unit: %
%% cutoff volume in PM2.5
    dV_cutoff_PM10=sum(dV(M_cutoff:M_PM10));
    percent_cutoff_PM10=dV_cutoff_PM10/dV_PM10*100; %unit: %


    
    
clc