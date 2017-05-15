% Copyright (C) 2013 Quan Wang <wangq10@rpi.edu>,
% Signal Analysis and Machine Perception Laboratory,
% Department of Electrical, Computer, and Systems Engineering,
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA

% dynamic time warping of two signals

function d=dLB(s,t,w)
% s: signal 1, size is ns*k, row for time, colume for channel 
% t: signal 2, size is nt*k, row for time, colume for channel 
% w: window parameter
%      if s(i) is matched with t(j) then |i-j|<=w
% d: resulting distance

if nargin<3
    w=Inf;
end

 ns=size(s,1);
 nt=size(t,1);
% if size(s,2)~=size(t,2)
%     error('Error in dtw(): the dimensions of the two input signals do not match.');
% end
% w=max(w, abs(ns-nt)); % adapt window size
%% LB_Keogh

U = zeros(ns,1);
L = zeros(ns,1);

for i = 1:ns
    U(i) = max(s(max(i-w,1):min(i+w,ns)));
    L(i) = min(s(max(i-w,1):min(i+w,ns)));
end

d = sqrt(sum([[t > U].* [t-U]; [t < L].* [L-t]].^2));
