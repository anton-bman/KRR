function allcoords = getGrid(gridRes, gridMinMax)
%gridRes = resolution of the grid will be gridRes*2+1
%gridMinMax = absolute value of min and max value of the grid
%dx = dy = dz = gridMinMax*2./(gridRes-1)

[x, y, z] = ndgrid(linspace(-gridMinMax, gridMinMax, gridRes), ... 
    linspace(-gridMinMax, gridMinMax, gridRes), ...
    linspace(-gridMinMax, gridMinMax, gridRes));

allcoords = x;
allcoords(:,:,:,2) = y;
allcoords(:,:,:,3) = z;

end