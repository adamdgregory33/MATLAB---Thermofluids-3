% Converts static to total temperature

function TTotal = TotalTempIdeal(T1, PT2, P1, y)
TTotal = T1.*(PT2./P1).^((y-1)/y);
end