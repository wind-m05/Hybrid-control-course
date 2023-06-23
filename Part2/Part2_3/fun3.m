function dxdt = fun3(t,x,Ap,Bp,Cp,Dp,a,r,d)
% state variable = [xp1 xp2 x1]'

xp1 = x(1);
xp2 = x(2);
x1 = x(3);

if -x1*Cp*[xp1; xp2] + x1*Dp*r >= 0
    dxdt = [Ap Bp; -Cp a]*x + [[0;0] Bp; Dp 0]*[r; d];
elseif -x1*Cp*[xp1; xp2] + x1*Dp*r < 0
    dxdt = [Ap -Bp; -Cp a]*x + [[0;0] Bp; Dp 0]*[r; d];
end

end