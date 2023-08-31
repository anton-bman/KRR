function [kernel_direc, kernel_nondirec, kernel_gauss] = kernelFuncGrid(... 
    grid, mic_coords, srcpos, k, beta, generalAngle)
%%%% Returns the kernel vector (matrix) as described in report
%%% interpCoords are the coordinates of the positions of where
%   we want to interpolate, in a grid format, as
%   interpCoords(:, 8, 24) would correspond to the
%   x-coordinates where y(yidx=8) and z(zidx=24)
%   --- Create a grid with ndgrid() or "help getGrid".
%%% mic_coords are the coordinates of the microphones of which
%   we have observed values, of size
%   [(1 x size(micx_coords,2)); 
%    (1 x size(micy_coords,2)); 
%    (1 x size(micz_coords,2))]
%%% srcpos are the coordinates of the sound (point) sources, of size
%   [(1 x size(src_xcoords,2)); 
%    (1 x size(src_ycoords,2)); 
%    (1 x size(src_ycoords,2))]
%%% k is the wavenumber and
%%% beta is the regularization parameter for the directional kernel
%%% generalAngle is set to true if you want the angle to the point source
%   to be calculated for each individual microphone and is set to false
%   if you want it to be based on the angle from the origin.


xinterpCoords = grid(:,:,:,1);
yinterpCoords = grid(:,:,:,2);
zinterpCoords = grid(:,:,:,3);



if not(generalAngle == true || generalAngle == false)
    generalAngle = true;
end


nmic = size(mic_coords,2);

for srcIdx = 1:size(srcpos,2)
    srcTemp = srcpos(:,srcIdx);
    if generalAngle == true
        [phi, thetas, ~] = cart2sph((srcTemp(1)), (srcTemp(2)), (srcTemp(3)));
        thetas = pi/2 - thetas;
    end
    for angleFrom = 1:nmic

        xr = xinterpCoords - mic_coords(1,angleFrom);
        yr = yinterpCoords - mic_coords(2,angleFrom);
        zr = zinterpCoords - mic_coords(3,angleFrom);
        absr = sqrt(xr.^2 + yr.^2 + zr.^2);
        
        if generalAngle == false
            print('need to fix the general angle for grid')
            [phi, thetas, ~] = cart2sph((srcTemp(1)-xinterpCoords(1,:)), (srcTemp(2)-yinterpCoords(2,:)), (srcTemp(3)-zinterpCoords(3,:)));
            thetas = pi/2 - thetas;
        end

        nu = [cos(phi).*sin(thetas); sin(phi).*sin(thetas);...
        cos(thetas)];
        nu_min_rvec(:,:,:,1) = xr;
        nu_min_rvec(:,:,:,2) = yr;
        nu_min_rvec(:,:,:,3) = zr;
        nu_min_rvec(:,:,:,1) = (1j*beta*nu(1) - k*nu_min_rvec(:,:,:,1)).^2;
        nu_min_rvec(:,:,:,2) = (1j*beta*nu(2) - k*nu_min_rvec(:,:,:,2)).^2;
        nu_min_rvec(:,:,:,3) = (1j*beta*nu(3) - k*nu_min_rvec(:,:,:,3)).^2; %.*conj(1j*beta*nu(3) - k*nu_min_rvec(:,:,:,3));
        
        if beta ~= 0
            kernel_direc(:,:,:,angleFrom,srcIdx) = 2*beta/(exp(beta)-exp(-beta))*sphBesselj(...
                sqrt(sum(nu_min_rvec,4)) );
        else
            kernel_direc(:,:,:,angleFrom,srcIdx) = sphBesselj(...
                sqrt(sum(nu_min_rvec,4)) );
        end
        kernel_nondirec(:,:,:,angleFrom) = sphBesselj(k*absr);
        %kernel_nondirec(:,:,:,angleFrom) = higherorder_kernel([xr, yr, zr], k);
        kernel_gauss(:,:,:,angleFrom) = exp(-(0.5*k)^2.*absr.^2);
    end
end