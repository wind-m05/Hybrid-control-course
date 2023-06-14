clear 
clc

%% Sliding mode lambda calculation
% syms lambda x1 x2
% eq = 4*x1-4*lambda*x1-5*x2+4*lambda*x2 == 2*x1*(-4*x1+6*lambda*x1-7*x2+13*lambda*x2);
% % eq = (2-2*lambda)-(5*x2/2*x1)+2*(x2/x1)*lambda == -4*x1+6*lambda*x1-7*x2+13*lambda*x2;
% eq_subs = subs(eq,x2,x1^2+1);
% lambda_sol_rody = solve(eq_subs,lambda)

%%
syms lambda x1 x2 dx1 dx2
dx1 = lambda*(2*x1 + 4*x2) + (1-lambda)*(-4*x1 - 7*x2);
dx2 = lambda*(0*x1 - 1*x2) + (1-lambda)*( 4*x1 - 5*x2) == 2*x1*dx1;

lambda_sol1 = solve(dx2,lambda);

dx2_sub = subs(dx2,x2,x1^2+1);
lambda_sol2 = solve(dx2_sub,lambda);