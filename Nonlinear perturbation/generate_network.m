function network = generate_network(R1,C1,S1,S2)
network = @(t,x,time,u) ODE(t,x,time,u,R1,C1,S1,S2);
end