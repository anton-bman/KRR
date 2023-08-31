function ROI_coords = ROI(allcoords, varargin)
%allcoords = grid coordinates (3d)
% ROI_coords = ROI(allcoords, varargin),
% where varagin has the following options:
% ROI_coords = ROI(allcoords, r, 'sphere')
% ROI_coords = ROI(allcoords, xSide, ySide, zSide, 'cube')
% where 'sphere' gives you a sphere with radius r, centered in the origin,
% and 'cube' gives you a cube with sides [xSide x ySide x zSide], also
% centered in the origin.




%extract indices
x =  allcoords(:,:,:,1);
y =  allcoords(:,:,:,2);
z =  allcoords(:,:,:,3);

string = varargin{end};

if string == 'sphere'
    r = varargin{1};
    ROI_coords = double(sqrt((x.^2 + y.^2 + z.^2)) < r);
elseif string == 'box'
    xSide = varargin{1};
    ySide = varargin{2};
    zSide = varargin{3};
    ROI_coords = double(abs(x) < xSide & abs(y) < ySide & abs(z) < zSide);
end

end







