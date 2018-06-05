function initConds = initCondsTNNP()

V = -84;
m = 0;
h = 1;
j = 1;
s = 1;
r = 0;
xs = 0;
xr1 = 0;
xr2 = 0;
d = 0;
f = 1;
f2 = 1;
fcass = 1;
Nai = 8;
Ki = 137;
Caitotal = 0.12e-3; 
Casrtotal = 0;
Casstotal = 0;
Rb = 1;

initConds = [V, m, h, j, s, r, xs, xr1, xr2, d, f, f2, fcass, Nai, ...
    Ki, Caitotal, Casrtotal, Casstotal, Rb];