function APD = findAPD(V, t, X)
% APD = findAPD(V, t, stimPeriod, X)
%
% Finds APDX

% First, find the derivative of the action potential
dV = diff(V)./diff(t);

% Now, find the end of the upstroke, which gives APA
APApos = find(dV < 0,1);
APA = V(APApos);

% Find what is the potential at X% repolarisation
VX = (1-X/100)*(APA-min(V)) + min(V);

% Now look for the point in V where that happens
APDX = find(V(APApos : end) <= VX , 1) + APApos - 1;

% Calculate APDX

APD = t(APDX) - t(APApos);
if isempty(APD)
    APD = 0;
end