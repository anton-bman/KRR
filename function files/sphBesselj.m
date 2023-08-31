function js = sphBesselj(x)
%%% Implementation of the 0th order spherical bessel function of
% the first kind.
% x can be in matrix form.

js = sinc(x/pi);

end