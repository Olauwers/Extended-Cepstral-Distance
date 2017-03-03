function d=dLBKeogh(s,t,w)

if nargin<3
    w=Inf;
end

 ns=size(s,1);
 nt=size(t,1);

%% LB_Keogh

U = zeros(ns,1);
L = zeros(ns,1);

for i = 1:ns
    U(i) = max(s(max(i-w,1):min(i+w,ns)));
    L(i) = min(s(max(i-w,1):min(i+w,ns)));
end

d = sqrt(sum([[t > U].* [t-U]; [t < L].* [L-t]].^2));
