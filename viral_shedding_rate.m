function y = viral_shedding_rate(t,p)
% viral shedding peaks at omega1/(2*omega2) when t=omega2

omega1     = p(1);
omega2     = p(2);

y = nan(length(t),1);

    y  = omega1*t./(omega2^2 + t.^2);

end