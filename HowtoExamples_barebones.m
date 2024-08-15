

%% Paramaters CsYbSe2 -- from our paper and extensive fitting. 
Bmn = [-0.4233    0.0117    0.5494    0.00035    0.0052   -0.00045];
Jxx = 0.5388; %Scaled From meV to K
Jzz = 0.6145;
X4 = [9.95 Bmn Jxx Jzz];    %
fac = [1,10,1,1e3,1e2,1e2]; % Scaling factor to scale steven's oerator energies equally for fitting. 

% This plots the CEF spectrum for the given Bmn Steven's operator
% coefficients. for H||ab and H||c
CEF_Spectrum_D3d('full spectrum',Bmn)

% Plots magnetic susceptibility for H||ab and H||c.
CEF_Magnetotropic_Modeling('Chi v T', X4);

%plots magnetotropic coefficient for H||ab and H||c. 
CEF_Magnetotropic_Modeling_Fast('principle k', [Bmn.*fac Jxx Jzz]);



% Bmn = [+0.4233    0.0117    0.5494    0.00035    0.2   -0.00045];
% Jxx = 0.5388; 
% Jzz = 0.6145;
% X4 = [9.95 Bmn Jxx Jzz];    %
% fac = [1,10,1,1e3,1e2,1e2]; % Scaling factor to scale steven's oerator energies equally for fitting. 
% CEF_Spectrum_D3d('full spectrum',Bmn)
% CEF_Magnetotropic_Modeling('Chi v T', X4);
% CEF_Magnetotropic_Modeling_Fast('principle k', [Bmn.*fac Jxx Jzz]);
