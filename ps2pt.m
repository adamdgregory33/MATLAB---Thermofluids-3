% calculates the total pressure from the static pressure, air properties
% and the velocity

function p_total = ps2pt(p_static, gamma, mach)

p_total = p_static*(1+((gamma-1)/2)*mach^2)^(gamma/(gamma-1));
end