% This function converts the static temperature to the total temperature

function t_total = ts2tt(ts, pt, p, gamma)
t_total = ts*(pt/p)^((gamma-1)/gamma);
end