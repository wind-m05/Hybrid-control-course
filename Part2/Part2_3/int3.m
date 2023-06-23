function dxdt = int3(t,x,Ap,Bp,Cp,Dp,a,r,d)
% state variable = [xp1 xp2 x1]'

dxdt = [Ap Bp; -Cp 0]*x + [[0;0] Bp; Dp 0] * [r; d];

end