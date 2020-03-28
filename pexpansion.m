% This function gives the pressure after expansion

function Pt5 = pexpansion(Pt4, Tt4, Tt5, gamma, efficiency)
Tdt5 = Tt4 - (Tt4-Tt5)/efficiency;

Pt5 = Pt4*(Tdt5/Tt4)^(gamma/(gamma-1));
end