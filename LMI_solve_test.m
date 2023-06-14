A = [-1, 2, 0;
     -3,-4, 1;
      0, 0,-2];
Q = 0.1*eye(3);
Pvar = sdpvar(3,3);
Lf = [A'*Pvar+Pvar*A+Q <= 0];
Lp = [Pvar >= 1e-9];
L = Lf + Lp;
diagnostics = optimize(L);
disp(diagnostics.info)
P = value(Pvar)

% Check values
eig(P)
A'*P+P*A+Q