% Add an input to the system.
function network_with_input = set_input(network,u,time)
network_with_input = @(t,x) network(t,x,time,u);
end