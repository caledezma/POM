function output=myTP06(time,X,flag_ode, param, cellType)

%% Give names to variables in X

V = X(1);
m = X(2);
h = X(3);
j = X(4);
s = X(5);
r = X(6);
xs = X(7);
xr1 = X(8);
xr2 = X(9);
d = X(10);
f = X(11);
f2 = X(12);
fcass = X(13);
Nai = X(14);
Ki = X(15);
Caitotal = X(16); 
Casrtotal = X(17);
Casstotal = X(18);
Rb = X(19);


%% Get external parameter values
Gna = param(1);
Gto = param(2);
Gks = param(3);
Gkr = param(4);
Gk1 = param(5);
Pnak = param(6);
Gcal = param(7);
knaca = param(8);
Gpca = param(9);
Gpk = param(10);
Gbna = param(11);
Gbca = param(12);

%% Define constants for the model
Cm = 0.185; %Cellular capacitance

R = 8314.472; %Gas constant
T = 310; %Temperature
F = 96485.3415; %Faraday constant

%Intracellular volumes
Vc = 0.016404; %Cytoplasmatic volume
Vsr = 0.001094; %Sarcoplasmatic reticulum volume
Vss = 0.00005468; %Subspace volume

%External ion concentrations
Nao = 140; %Extracellular Na+ concentration
Ko = 5.4; %Extracellular K+ concentration;
Cao = 2; %Extracellular Ca2+ concentration

%Calcium buffering dynamics constants
Kbufc = 0.001; %Cai half saturation constant for cytoplasmatic buffer
Bufc = 0.2; %Total cytoplasmatic buffer concentration
Kbufsr = 0.3; %Casr half saturation constant for sarcoplasmatic buffer
Bufsr = 10; %Total sarcoplasmatic buffer concentration
Kbufss = 0.00025; %Cass half saturation constant for subspace buffer
Bufss = 0.4; %Total subspace buffer concentration
pkna = 0.03; %Relative Iks permeability to Na+

% Ionic current constants
Kmk = 1; %Ko half saturarion constant of Inak
Kmna = 40; %Nai half saturation constant of Inak
gamma = 0.35; %Voltage dependence parameter of Inaca
alpha = 2.5; %Factor enhancing outward nature of Inaca
Kmnai = 87.5; %Nai half saturation constant for Inaca
Kmca = 1.38; %Cai half saturation constant for Inaca
ksat = 0.1; %Saturation factor for Inaca
Kpca = 0.0005; %Half saturation constant for Ipca

%Intracellular calcium flux dynamics parameters
Vleak = 0.00036; %Maximal Ileak conductance
Vmaxup = 0.006375; %Maximal Iup consuctance
Kup = 0.00025; %Half saturation constant for Iup
Vrel = 0.102; %Maximal Irel conductance
Vxfer = 0.0038; %Maximal Ixfer conductance
k1p = 0.15; %R to O and RI to I Irel transition rate
k2p = 0.045; %O to I and R to RI Irel transition rate
k3 = 0.06; %O to R and I to RI Irel transition rate
k4 = 0.005; %I to O and RI to I Irel transition rate
maxsr = 2.5; %Maximum value of kcasr
minsr = 1; %Minimum value of kcasr
EC = 1.5; %Casr half saturation constant of kcasr



%% Solve the equations

% Calcium transient values

Cai = 1/2 * (-(Kbufc+Bufc-Caitotal)+sqrt((Kbufc+Bufc-Caitotal).^2 + 4*Caitotal*Kbufc));
Casr = 1/2 * (-(Kbufsr+Bufsr-Casrtotal)+sqrt((Kbufsr+Bufsr-Casrtotal).^2 + 4*Casrtotal*Kbufsr));
Cass = 1/2 * (-(Kbufss+Bufss-Casstotal)+sqrt((Kbufss+Bufss-Casstotal).^2 + 4*Casstotal*Kbufss));

%Calculate reverse potentials
Ena = R*T / F * log(Nao*Nai.^-1);
Ek = R*T / F * log(Ko*Ki.^-1);
Eks = R*T / F * log((Ko + pkna*Nao)*(Ki + pkna*Nai).^-1);
Eca = R*T / (2*F) * log(Cao*Cai.^-1);

%Calculate Ion currents
% Fast Na+ current
Ina = Gna .* m.^3 .* h .* j .* (V - Ena);

% Transient outward current
Ito = Gto .* r .* s .* (V - Ek);

% Slow delayed rectifier current
Iks = Gks .* xs.^2 .* (V - Eks);

% Rapid delayed rectifier current
Ikr = Gkr * sqrt(Ko/5.4) .* xr1 .* xr2 .* (V - Ek);

% Inward rectifier K+ current
alphak1 = 0.1 * (1 + exp(0.06*(V-Ek-200))).^-1;
betak1 = (3*exp(0.0002*(V-Ek+100)) + exp(0.1*(V-Ek-10))) .* ...
    (1 + exp(-0.5*(V-Ek))).^-1;
xk1inf = alphak1 ./ (alphak1 + betak1);
Ik1 = Gk1 * xk1inf .* (V - Ek) ;

% Na+/K+ pump current
Inak = Pnak * (Ko * Nai) .* ((Ko+Kmk) * (Nai+Kmna) .* ...
    (1 + 0.1245*exp(-0.1*V*F/(R*T)) + 0.0353*exp(-V*F/(R*T)))).^-1;

% L-type Ca2+ current
Ical = Gcal * d .* f .* f2 .* fcass * 4 .* (V-15)*F^2/(R*T) .* ...
    (0.25*Cass.*exp(2*(V-15)*F/(R*T))-Cao) .* (exp(2*(V-15)*F/(R*T))-1).^-1;

% Na+/Ca2+ exchanger current
Inaca = knaca * (exp(gamma*V*F/(R*T)).*Nai.^3*Cao - ...
    exp((gamma-1)*V*F/(R*T)).*Nao^3.*Cai*alpha) .* ((Kmnai^3+Nao^3) * ...
    (Kmca+Cao) .* (1+ksat*exp((gamma-1)*V*F/(R*T)))).^-1;

% Plateau currents
Ipca = Gpca * Cai .* (Kpca + Cai).^-1;
Ipk = Gpk * (V - Ek) .* (1 + exp((25-V)/5.98)).^-1;

% Background currents
Ibna = Gbna * (V - Ena);
Ibca = Gbca * (V - Eca);

%Calculate intracellular calcium dynamics
kcasr = maxsr - ((maxsr - minsr) .* (1+(EC./Casr).^2).^-1);
k1 = k1p*kcasr.^-1;
k2 = k2p*kcasr;
O = k1 .* Cass.^2 .* Rb ./ (k3 + k1.*Cass.^2);

Ileak = Vleak .* (Casr - Cai);
Iup = Vmaxup .* (1 + Kup^2./Cai.^2).^-1;
Irel = Vrel * O .* (Casr-Cass);
Ixfer = Vxfer * (Cass - Cai);

%Calculate variations on gating variables

% Ina gates
minf = (1 + exp((-56.86 - V)/9.03)).^-2;
alpham = (1 + exp((-60 - V)/5)).^-1;
betam = 0.1 * ((1 + exp((V+35)/5)).^-1 + (1 + exp((V - 50)/200)).^-1);
taum = alpham .* betam;

hinf = (1 + exp((V + 71.55)/7.43)).^-2;
if V >= -40
    alphah = 0;
    betah = 0.77 * (0.13 * (1 + exp(-(V+10.66)/11.1))).^-1;
else
    alphah = 0.057*exp(-(V+80)/6.8);
    betah = 2.7 * exp(0.079*V) + 3.1e5 * exp(0.3485*V);
end
tauh = (alphah + betah).^-1;

jinf = (1 + exp((V+71.55)/7.43)).^-2;
if V >= -40
    alphaj = 0;
    betaj = 0.6 .* exp(0.057*V) .* (1 + exp(-0.1*(V+32))).^-1;
else
    alphaj = (-2.5428e4*exp(0.2444*V)-6.948e-6*exp(-0.04391*V)).*(V + 37.78)...
        .*(1 + exp(0.311*(V + 79.23))).^-1;
    betaj = 0.02424*exp(-0.01052*V).*(1+exp(-0.1378*(V+40.14))).^-1;
end
tauj = (alphaj + betaj).^-1;

% Ito gates
rinf = (1 + exp((20-V)/6)).^-1;
taur = 9.5 * exp(-(V+40).^2/1800) + 0.8;

if cellType == 1
    sinf = (1 + exp((V+28)/5)).^-1;
    taus = 1000 * exp(-(V+67).^2/1000) + 8;
else
    sinf = (1 + exp((V+20)/5)).^-1;
    taus = 85 * exp(-(V+45).^2/320) + 5 * (1 + exp((V-20)/5)).^-1 + 3;
end

% Iks gates
xsinf = (1 + exp((- 5 - V)/14)).^-1;
alphaxs = 1400 * (sqrt(1 + exp((5 - V)/6))).^-1;
betaxs = (1 + exp((V - 35)/15)).^-1;
tauxs = alphaxs .* betaxs + 80;

% Ikr gates
xr1inf = (1 + exp((- 26 - V)/7)).^-1;
alphaxr1 = 450 * (1 + exp((- 45 - V)/10)).^-1;
betaxr1 = 6 * (1 + exp((V + 30)/11.5)).^-1;
tauxr1 = alphaxr1 .* betaxr1;

xr2inf = (1 + exp((V + 88)/24)).^-1;
alphaxr2 = 3 * (1 + exp((- 60 - V)/20)).^-1;
betaxr2 = 1.12 * (1 + exp((V - 60)/20)).^-1;
tauxr2 = alphaxr2 .* betaxr2;

% Ical gates
dinf = (1 + exp((- 8 - V)/7.5)).^-1;
alphad = 1.4 * (1 + exp((- 35 - V)/13)).^-1 + 0.25;
betad = 1.4 * (1 + exp((V + 5)/5)).^-1;
gammad = (1 + exp((50-V)/20)).^-1;
taud = alphad .* betad + gammad;

finf = (1 + exp((V + 20)/7)).^-1;
alphaf = 1102.5 * exp(-((V+27)/15).^2);
betaf = 200 * (1 + exp((13 - V)/10)).^-1;
gammaf = 180 * (1 + exp((V+30)/10)).^-1 + 20;
tauf = alphaf + betaf + gammaf;

f2inf = 0.67 * (1 + exp((V + 35)/7)).^-1 + 0.33;
alphaf2 = 600 * exp(-(V + 25).^2/170);
betaf2 = 31 * (1 + exp((25-V)/10)).^-1;
gammaf2 = 16 * (1 + exp((V + 30)/10)).^-1;
tauf2 = alphaf2 + betaf2 + gammaf2;

fcassinf = 0.6 * (1 + (Cass/0.05).^2).^-1 + 0.4;
taufcass = 80 * (1 + (Cass/0.05).^2).^-1 + 2;

% Find variations in gates
dm = (minf - m) ./ taum;
dh = (hinf - h) ./ tauh;
dj = (jinf - j) ./ tauj;
dr = (rinf - r) ./ taur;
ds = (sinf - s) ./ taus;
dxs = (xsinf - xs) ./ tauxs;
dxr1 = (xr1inf - xr1) ./ tauxr1;
dxr2 = (xr2inf - xr2) ./ tauxr2;
dd = (dinf - d) ./ taud;
df = (finf - f) ./ tauf;
df2 = (f2inf - f2) ./ tauf2;
dfcass = (fcassinf - fcass) ./ taufcass;

% Calculate stimulation current
amp = -80.0;
dur = 0.5;
if time <= dur
    Istim = amp;
else
    Istim = 0;
end

%Calculate variations in intracellular ion concentrations
%Sodium
dNai = -(Ina + Ibna + 3*Inak + 3*Inaca)/(Vc*F)*Cm;

%Potasium
dKi = - (Ik1 + Ito + Ikr + Iks - 2*Inak + Ipk + Istim) * Cm / (Vc * F);

%Calcium
dCaitotal = -(Ibca+Ipca-2*Inaca) * Cm / (2*Vc*F) + (Vsr/Vc) * (Ileak-Iup) + Ixfer;
dCasrtotal = Iup - Ileak - Irel;
dCasstotal = -Ical*Cm/(2*Vss*F) + Vsr*Irel/Vss - Vc*Ixfer/Vss;

dRb = k4*(1-Rb) - k2.*Cass.*Rb;

% Variation in voltage
dV = -(Ina+ Ito+Iks+Ikr+Ik1+Inak+Ical+Inaca+Ipca+Ipk+Ibna+Ibca+Istim);

%% Return results
if flag_ode==1
    output=[dV, dm, dh, dj, ds, dr, dxs, dxr1, dxr2, dd, df, df2, dfcass,...
        dNai, dKi, dCaitotal, dCasrtotal, dCasstotal, dRb]';
else
    output=[Ina, Ito, Iks, Ikr, Ik1, Inak, Ical, Inaca, Ipca, Ipk, Ibna, Ibca];
end

