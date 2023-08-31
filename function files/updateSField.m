function z_upd = updateSField(A,r, k, w, t)
%Calculates the sound pressure at distance r from sound source

z_upd = (A./(r)).*(exp(-1i*(k*(r) - w*t)));
%z_upd = A*exp(-1i*(k*(r) - w*t));
end