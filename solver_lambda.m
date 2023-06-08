clear all, close all, clc

%% Sliding mode lambda calculation
syms lambda x1 x2
eq = 4*x1-4*lambda*x1-5*x2+4*lambda*x2 == 2*x1*(-4*x1+6*lambda*x1-7*x2+13*lambda*x2);
% eq = (2-2*lambda)-(5*x2/2*x1)+2*(x2/x1)*lambda == -4*x1+6*lambda*x1-7*x2+13*lambda*x2;
eq_subs = subs(eq,x2,x1^2+1);
lambda_sol_rody = solve(eq_subs,lambda)
