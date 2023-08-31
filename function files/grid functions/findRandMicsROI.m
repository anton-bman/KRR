function [mic_coords, mic_idx] = findRandMicsROI(allcoords, nmic, ROI_coord)
%allcoords = grid coordinates (3d)
%r = radius of the region of interest (ROI)
%nmic = number of microphones in the ROI
%ROI_coords = coordinates of the ROI, (x-coord; y-coord; z-coord)
%%% Returns microphones chosen uniformly in the ROI.

%extract indices
x =  allcoords(:,:,:,1);
y =  allcoords(:,:,:,2);
z =  allcoords(:,:,:,3);

%extract the values
sz = size(x);
x_calc =  x(:,1,1);
y_calc =  y(1,:,1)';
z_calc =  reshape(z(1,1,:), [1,sz(3)]);

%find the indices
[I,J,K] = ind2sub(size(ROI_coord),find(ROI_coord==1));

%pick nmic random indices
mic_idx = randperm(length(I),nmic);

%return to values instead of indices
mic_coords = [x_calc(I(mic_idx))'; y_calc(J(mic_idx))' ; z_calc(K(mic_idx)')];


end