function L = SL(i,max,a,b)
%SL This function simulates a non-linear saturating inductor.
%   i represents the current that flows through the inductor, max is the
%   maximum value of inductance L, b is the value at which the inductor
%   starts saturating, a is the value at which the inductor is completely
%   saturated.
L = 1;
if (i > -a) && (i <= -b)
    L = 1+((max-1)/abs(a-b))*(a+i);
elseif (i > -b) && (i <= b)
    L = max;
elseif (i > b) && (i <= a)
    L = 1 + ((max-1)/abs(a-b))*(a-i);
end
end
