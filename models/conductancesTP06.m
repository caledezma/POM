function params = TNNPparams(cellType)
% cellType = 1 (endo), cellType = 2 (M-cells), cellType = 3 (epi) 

Gna = 14.838; %Maximum Ina conductance
%Maximal Ito conductance
if cellType == 1
    Gto = 0.073; %Endocardial
else
    Gto = 0.294; %Epicardial
end
%Maximal Iks conductance
if cellType == 2 
    Gks = 0.098; %For M cells
else
    Gks = 0.392; %For epi and endocardial cells
end
Gkr = 0.153; %Maximal Ikr conductance
Gk1 = 5.405; %Maximal Ik1 conductance
Pnak = 2.724; %Maximal Inak
Gcal = 3.980e-5; %Maximum Ical conductance
knaca = 1000; %Maximal Inaca
Gpca = 0.1238; % Maximum Ipca conductance
Gpk = 0.0146; %Maximum Ipk conductance
Gbna = 0.00029; %Maximum Ibna conductance
Gbca = 0.000592; %Maximum Ibca conductance

params = [Gna , Gto , Gks , Gkr , Gk1 , Pnak , Gcal , ...
    knaca , Gpca , Gpk , Gbna , Gbca];