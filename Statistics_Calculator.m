clc; clear all;

%% %%%%%%%%%%%%%%%%%%%%% Ultimate Statistics f(x,y) Calculator %%%%%%%%%%%%%%%%%%%%%%%%

syms f(x,y) g(x,y)
syms a b x y k n t
f(x,y) = 2-x-y;
g(x,y) = a*x + b*y;
h(t) = exp(t) * (1-exp(n*t)) / (n*(1-exp(t)));

Xbounds = [0 1]; % more complex inequality
Ybounds = [0 1];
X_var_bounds = [0 1]; % bounds for individual variables
Y_var_bounds = [0 1];
givenX = 0.5;
givenY = 0.5;
a = 1;
b = 1;
testBounds = [0 0.25];
%testBounds = [0.8 1];
% diff(f, x, t) is derivative of f with respect to x at t
% int(f, x, bounds) is integral of f with respect to x with [bounds]

Fx = int(f, y, Ybounds);
Fy = int(f, x, Xbounds);
Fx(x,y) = Fx;
Fy(x,y) = Fy;
PFxbounds = int(Fx, x, testBounds);
PFybounds = int(Fy, y, testBounds);
Fyx = simplify(f / Fx);
Fxy = simplify(f / Fy);
muYXgiven = int(Fyx * y, y, Ybounds);
muXYgiven = int(Fxy * x, x, Xbounds);
muY2Xgiven = int(Fyx * y^2, y, Ybounds); % only used for LOTUS calc of variance of conditionals with given
muX2Ygiven = int(Fxy * x^2, x, Xbounds);
varYXgiven = muY2Xgiven - muYXgiven^2;
varXYgiven = muX2Ygiven - muXYgiven^2;
PYboundsXgivenFun = int(Fyx, y, testBounds);
PXboundsYgivenFun = int(Fxy, x, testBounds);

% resets multivariable to single variable functions
Fx = int(f, y, Ybounds);
Fy = int(f, x, Xbounds);
muX = int((x * Fx), x, X_var_bounds);
muY = int((y * Fy), y, Y_var_bounds);
muX2 = int((x^2 * Fx), x, X_var_bounds);
muY2 = int((y^2 * Fy), y, Y_var_bounds);
muXY = int(int(x*y*f, y, Ybounds), x, X_var_bounds);
varX = muX2 - (muX)^2;
varY = muY2 - (muY)^2;
stdX = sqrt(varX);
stdY = sqrt(varY);
covXY = muXY - muX*muY;
corrXY = covXY / (sqrt(varX * varY));
muG = a*muX + b*muY;
var_aXbY = a^2 * varX + b^2 * varY + 2*covXY;

fprintf("F(x) = " + string(Fx) + "\n");
fprintf("F(y) = " + string(Fy) + "\n");
fprintf("P(" + testBounds(1) + " < X < " + testBounds(2) + ") = " + string(PFxbounds) + " = " + double(PFxbounds) +  "\n");
fprintf("P(" + testBounds(1) + " < Y < " + testBounds(2) + ") = " + string(PFybounds) + " = " + double(PFybounds) +  "\n");

fprintf("\n");
fprintf("mu(x) = " + string(muX) + " = " + double(muX) +  "\n");
fprintf("mu(y) = " + string(muY) + " = " + double(muY) +  "\n");
fprintf("mu(x^2) = " + string(muX2) + " = " + double(muX2) +  "\n");
fprintf("mu(y^2) = " + string(muY2) + " = " + double(muY2) +  "\n");
fprintf("mu(xy) = " + string(muXY) + " = " + double(muXY) +  "\n");
fprintf("\n");
fprintf("var(x) = " + string(varX) + " = " + double(varX) +  "\n");
fprintf("var(y) = " + string(varY) + " = " + double(varY) +  "\n");
fprintf("std(x) = " + string(stdX) + " = " + double(stdX) +  "\n");
fprintf("std(y) = " + string(stdY) + " = " + double(stdY) +  "\n");
fprintf("\n");
fprintf("cov(x,y) = " + string(covXY) + " = " + double(covXY) + "\n");
fprintf("corr(x,y) = " + string(corrXY) + " = " + double(corrXY) + "\n");
fprintf("\n");
fprintf("f(Y|X) = " + string(Fyx) + "\n");
fprintf("f(X|Y) = " + string(Fxy) + "\n");
fprintf("E(Y|X) = " + string(muYXgiven) + "\n");
fprintf("E(X|Y) = " + string(muXYgiven) + "\n");
fprintf("E(Y|X=" + givenX + ") = " + string(muYXgiven(givenX)) + " = " + double(muYXgiven(givenX)) + "\n");
fprintf("E(X|Y=" + givenY + ") = " + string(muXYgiven(givenY)) + " = " + double(muXYgiven(givenY)) + "\n");
fprintf("var(X|Y) = " + string(varXYgiven) + "\n");
fprintf("var(Y|X) = " + string(varYXgiven) + "\n");
fprintf("var(X|Y=" + givenY + ") = " + string(varXYgiven(givenY)) + "\n");
fprintf("var(Y|X=" + givenX + ") = " + string(varYXgiven(givenX)) + "\n");

fprintf("\n");
fprintf("P(" + testBounds(1) + " < Y < " + testBounds(2) + "|X) = " + string(PYboundsXgivenFun) + "\n");
fprintf("P(" + testBounds(1) + " < X < " + testBounds(2) + "|Y) = " + string(PXboundsYgivenFun) + "\n");
fprintf("P(" + testBounds(1) + " < Y < " + testBounds(2) + "|X = " + givenX + ") = " + string(PYboundsXgivenFun(givenX)) + " = " + double(PYboundsXgivenFun(givenX)) + "\n");
fprintf("P(" + testBounds(1) + " < X < " + testBounds(2) + "|Y = " + givenY + ") = " + string(PXboundsYgivenFun(givenY)) + " = " + double(PXboundsYgivenFun(givenY)) + "\n");

fprintf("\n");
fprintf("E[g(x,y)] = " + string(muG) + " = " + double(muG) + "\n");
fprintf("Var[g(x,y)] = " + string(var_aXbY) + " = " + double(var_aXbY) + "\n");

