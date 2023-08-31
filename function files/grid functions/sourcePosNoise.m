function srcpos_noise = sourcePosNoise(srcpos, noiseError)
% noiseError is the expected value of the error in DEGREES

degreeError = noiseError*pi/180;
for idx = 1:size(srcpos,2)
    srcTemp = srcpos(:,1);
    [phi, thetas, rs] = cart2sph((srcTemp(1)), (srcTemp(2)), (srcTemp(3)));
    thetas = thetas + normrnd(0,degreeError*sqrt(pi/2));
    phi = phi + normrnd(0,degreeError*sqrt(pi/2));

    [srcpos_noise(1,idx), srcpos_noise(2,idx), srcpos_noise(3,idx)] = sph2cart(phi, thetas, rs);
end

end