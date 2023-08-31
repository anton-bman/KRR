function [kernel_direc, kernel_nondirec, kernel_gauss] = kernelFunc(...
    interpCoords, mic_coords, srcpos, k, beta, generalAngle)
%%% If using a grid for interpCoords, "help kernelFuncGrid"
%%% If using no grid for interpCoords, "help kernelFuncNoGrid".

if length(size(interpCoords)) == 4

    [kernel_direc, kernel_nondirec, kernel_gauss] = kernelFuncGrid(...
        interpCoords, mic_coords, srcpos, k, beta, generalAngle);

elseif length(size(interpCoords)) == 2
    [kernel_direc, kernel_nondirec, kernel_gauss] = kernelFuncNoGrid(...
        interpCoords, mic_coords, srcpos, k, beta, generalAngle);
end

end