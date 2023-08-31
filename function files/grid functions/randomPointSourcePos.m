function srcpos = randomPointSourcePos(r, r2,  ...
    nSrc)
% r = Minimum radius from the center of the grid
% r2 = Maximum radius added to the r, which leads
%   to the sources being between [r,r+r2] from the
%   center of the grid
% nSrc = Number of point sources emitting sound waves

rho = rand(nSrc,1)*r2+r;
thetas = rand(nSrc,1)*pi;
phi = rand(nSrc,1)*2*pi;
xpos = rho.*sin(thetas).*cos(phi);
ypos = rho.*sin(thetas).*sin(phi);
zpos = rho.*cos(thetas);
srcpos = [xpos';ypos';zpos'];
end