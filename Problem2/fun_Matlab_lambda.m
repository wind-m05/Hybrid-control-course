function tmp = fun_Matlab_lambda(t,x)

x1 = x(1);
x2 = x(2);

A1 = [2 4; 0 -1];
A2 = [-4 -7; 4 -5];

eps = 1e-3;

if x2 > x1^2 + 1 + eps
    tmp = A1*x;
elseif x2 < x1^2 + 1 - eps
    tmp = A2*x;
else %if x2 == x1^2 + 1
    lambda = (14*x1^3+3*x1^2+18*x1-5)/(22*x1^3+8*x1^2+26*x1-4);
    tmp = lambda*A1*x + (1-lambda)*A2*x;
end