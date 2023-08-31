function [pHhzDir, pHhz, pGauss] = interpAllGrid(...
    pressure, nmic, lambda, ...
    KDir, KNoDir, KGauss, kappaDir, kappaNonDir, kappaGauss)
%%%
%[pHhzDir, pHhz, pGauss] = pressure_interp_real(...
%    pressure, nmic, lambda, ...
%    KDir, KNoDir, KGauss, kappaDir, kappaNonDir, kappaGauss)
%
%%% Interpolates the sound field for
% 1. the directional helmholtz kernel (pHhzDir)
% 2. the non-directional helmholtz kernel (pHhz)
% 3. the Gaussian kernel (pGauss)
%
% pressure = the pressure recorded in the microphones
% nmic = number of microphones
% lambda = regularization parameter
% K = K-matrix as defined in paper (Kdir, KNoDir, KGauss)
% kappa = kernel as defined in paper (kappaDir, kappaNonDir, kappaGauss)
%
%%% The K-matrix must be 2-dimensional, and the kappa-matrix must be
% 4-dimensional for the grid, and 2-dimensional for the non-grid.

%%%
nx = size(kappaDir,1);
ny = size(kappaDir,2);
nz = size(kappaDir,3);

%the matrix ((K+lambda*I_M)^(-1) s)^T for the new kernel
kMatDir = (inv(KDir + lambda.*eye(size(KDir))) * pressure.').';

%the matrix ((K+lambda*I_M)^(-1) s)^T for the gaussian kernel
kMatGauss = (inv(KGauss + lambda.*eye(size(KGauss))) * pressure.').';

%the matrix ((K+lambda*I_M)^(-1) s)^T for the nondir kernel
kMatNonDir = (inv(KNoDir + lambda.*eye(size(KNoDir))) * pressure.').';

if length(size(kappaDir)) == 4 %if we are using a grid

    %%%this is the part that takes the most time -- could be reduced by
    %%%changing the kappa matrix when implementing it instead
    %directional kernel
    D = permute(kappaDir, [4 1 2 3]);
    E = reshape(D,nmic, nx*ny*nz);
    F = kMatDir*E;

    pHhzDir = reshape(F,nx,ny,nz);

    %gaussian kernel
    D_gauss = permute(kappaGauss, [4 1 2 3]);
    E_gauss = reshape(D_gauss,nmic, nx*ny*nz);
    F_gauss = kMatGauss*E_gauss;

    pGauss = reshape(F_gauss,nx,ny,nz);


    %nondirectional kernel
    D_old = permute(kappaNonDir, [4 1 2 3]);
    E_old = reshape(D_old,nmic, nx*ny*nz);
    F_old = kMatNonDir*E_old;

    pHhz = reshape(F_old,nx,ny,nz);

elseif length(size(kappaDir)) == 2
    %this is for the case where we are using the coordinates,
    % and not a grid
    pHhzDir = kMatDir*kappaDir.';
    pHhz    = kMatNonDir*kappaNonDir.';
    pGauss  = kMatGauss*kappaGauss.';
end

end