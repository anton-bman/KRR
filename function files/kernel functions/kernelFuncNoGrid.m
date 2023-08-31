function [kernel_direc, kernel_nondirec, kernel_gauss] = kernelFuncNoGrid(... 
    interpCoords, mic_coords, srcpos, k, beta, generalAngle)
%%%% Returns the kernel vector (matrix) as described in report
%%% interpCoords are the coordinates of the microphones of which
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

meanpoint = mean(mic_coords,2);

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
        xr = interpCoords(1,:) - mic_coords(1,angleFrom);
        yr = interpCoords(2,:) - mic_coords(2,angleFrom);
        zr = interpCoords(3,:) - mic_coords(3,angleFrom);
        absr = sqrt(xr.^2 + yr.^2 + zr.^2);
        
        if generalAngle == false
            [phi, thetas, ~] = cart2sph((srcTemp(1)-interpCoords(1,1)), (srcTemp(2)-interpCoords(2,1)), (srcTemp(3)-interpCoords(3,1)));
            thetas = pi/2 - thetas;
        end

        nu = [cos(phi)*sin(thetas); sin(phi)*sin(thetas);...
        cos(thetas)];
        nu_min_rvec(:,1) = xr;
        nu_min_rvec(:,2) = yr;
        nu_min_rvec(:,3) = zr;
        nu_min_rvec(:,1) = (1j*beta*nu(1) - k*nu_min_rvec(:,1)).^2;
        nu_min_rvec(:,2) = (1j*beta*nu(2) - k*nu_min_rvec(:,2)).^2;
        nu_min_rvec(:,3) = (1j*beta*nu(3) - k*nu_min_rvec(:,3)).^2; %.*conj(1j*beta*nu(3) - k*nu_min_rvec(:,:,:,3));
        
        if beta ~= 0
            kernel_direc(:,angleFrom,srcIdx) = 2*beta/(exp(beta)-exp(-beta))*sphBesselj(...
                sqrt(sum(nu_min_rvec.')) );
        else
            kernel_direc(:,angleFrom,srcIdx) = sphBesselj(...
                sqrt(sum(nu_min_rvec.')) );
        end
        kernel_nondirec(:,angleFrom) = sphBesselj(k*absr);
        %kernel_nondirec(:,angleFrom) = higherorder_kernel([xr, yr, zr], k);
        kernel_gauss(:,angleFrom) = exp(-(0.5*k)^2.*absr.^2);
    end
end