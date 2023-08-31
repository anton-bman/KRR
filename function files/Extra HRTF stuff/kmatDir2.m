function [K, K_nodir, K_gauss] = kmatDir2(mic_coords, k, srcpos, beta, generalAngle)
%returns the K-matrix as defined in the paper
% mic_coords = coordinates for the microphones
% k = wavenumber
% general angle... just ignore it, set to true and be happy

nmic = size(mic_coords,2);

K = zeros(nmic);
K_nodir = zeros(nmic);
K_gauss = zeros(nmic);

for srcposIdx = 1:size(srcpos,2)
    srcTemp = srcpos(:,srcposIdx);
    [phi, thetas, ~] = cart2sph((srcTemp(1)), (srcTemp(2)), (srcTemp(3)));
    thetas = pi/2 - thetas;
    if generalAngle == true
        [phi, thetas, ~] = cart2sph((srcTemp(1)), (srcTemp(2)), (srcTemp(3)));
        thetas = pi/2 - thetas;
    end
    for fromMic = 1:nmic
        %cartesian to spherical coordinate stuff
        %         [phi, thetas, ~] = cart2sph((srcTemp(1) - mic_coords(1,fromMic)), (srcTemp(2) - mic_coords(2,fromMic)), (srcTemp(3) - mic_coords(1,fromMic)));
        %         thetas = pi/2 - thetas;
        if generalAngle == false
            [phi, thetas, ~] = cart2sph((srcTemp(1)-mic_coords(1,fromMic)), (srcTemp(2)-mic_coords(1,fromMic)), (srcTemp(3)-mic_coords(1,fromMic)));
            thetas = pi/2 - thetas;
        end

        for toMic = 1:nmic

            xr = (mic_coords(1,fromMic)-mic_coords(1,toMic));
            yr = (mic_coords(2,fromMic)-mic_coords(2,toMic));
            zr = (mic_coords(3,fromMic)-mic_coords(3,toMic));
            absr = sqrt(xr^2 + yr^2 + zr^2);

            nu = [cos(phi)*sin(thetas); sin(phi)*sin(thetas);...
                cos(thetas)];
            rvec = [xr; yr; zr];
            if beta ~= 0
                K(fromMic,toMic,srcposIdx) = 2*beta/(exp(beta)-exp(-beta))*sphBesselj(...
                    sqrt( sum((1j.*beta*nu - k*rvec).^2) ));
            else
                K(fromMic,toMic,srcposIdx) = sphBesselj(...
                    sqrt( sum((1j.*beta*nu - k*rvec).^2) ));
            end
            K_nodir(fromMic,toMic) = sphBesselj(k.*absr);
            %K_nodir(fromMic,toMic) = higherorder_kernel([xr, yr, zr], k);
            K_gauss(fromMic,toMic) = exp(-(0.5*k).^2.*absr.^2);
        end
    end
end

end
