% This function gives the temperature after compression

function Tt2 = tcompression(Tt1, Pt1, Pt2, gamma, efficiency)

Tdt2 = Tt1*(Pt2/Pt1)^((gamma-1)/gamma);

Tt2 = Tt1 + (Tdt2-Tt1)/efficiency;
end