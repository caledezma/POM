function isPhysio = isPhysio(V,t,APD)
% Find out if an AP is physiological

% Derivative of V
dV = diff(V)./diff(t);

% Upstroke duration
UPD = t(find(dV < 0, 1));

if APD > 100 && APD < 400 && max(V) > 0 && min(V) < -64 && UPD < 10
    isPhysio = true;        
else
    isPhysio = false;
end
