% Dynamical equations of the system.
function dX = ODE(t,x,time,u,R1,C1,S1,S2)
uint = interp1(time,u,t);
dX = zeros(4,1);
dX(1) = (1/S1(x(1)))*x(2);
dX(2) = -(1/C1)*x(1)-(1/C1)*x(3) + (1/C1)*uint;
dX(3) = (1/S2(x(3)))*x(2) - (R1/S2(x(3)))*x(3);
dX(4) = dX(2) - R1*dX(3);
end