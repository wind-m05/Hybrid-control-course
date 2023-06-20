function tmp = fun_Matlab(t,x)

x1 = x(1);
x2 = x(2);

A1 = [2 4; 0 -1];
A2 = [-4 -7; 4 -5];

if x2 >= x1^2 + 1
    tmp = A1*x;
elseif x2 < x1^2 + 1
    tmp = A2*x;
end