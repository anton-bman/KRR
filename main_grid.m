clear
clc


%%%~~~~~~~------- defining the grid -------~~~~~~~%%%
gridRes = 41; %     <---- grid_res^2+1 is the res of grid
gridMinMax = 1;%    <---- min and max value of grid
allcoords = getGrid(gridRes, gridMinMax); %<- the grid

%distance between the points on the grid (in one direction)
dx = gridMinMax*2./(gridRes-1);
%%%~~~~~~~---- finding ROI coordinates ----~~~~~~~%%%
r = 0.5; %radius of the ROI, if spherical ROI
%  xSide = 1; this is an alternative, cuboid shape
%  ySide = 1;
%  zSide = 1;

ROI_coords = ROI(allcoords, r, 'sphere');
% ROI_coord = ROI_box(allcoords, xSide, ySide, zSide, 'box');

%%%~~~~~~~------ frequencies used etc -----~~~~~~~%%%
freq = [150, 200, 250, 300, 350, 400, 450, 500, 700, 900, 1100, 1500, 2000];%, 2500, 3000, 3500, 4000];
sound_speed = 343; %roughly 343 at 20 degree C

%%%~~~~~~~------ variables to define ------~~~~~~~%%%
beta = 6; %regularization parameter for the direction
lambda = 10^(-2); %regularization parameter
mic_stdev = 0;%0.04*sqrt(pi/2)/2; % 0.00*sqrt(pi/2)/2; %expected value of the error in m
noiseError = 0; %how many degrees of an error (expected value) there is in the estimated angles

%dB scale, SNR (signal-to-noise ratio)
%set to -inf to remove it completely
pressureNoise = -inf;

nIter = 4; %number of iterations
nmic = 12;
nSrc = 1; %number of point sources emitting sound

%%%~~~~~~~-------- wave properties --------~~~~~~~%%%
A = ones(1,nSrc);                  %Source strength
w = 2*pi.*freq;          %angular frequency
lambda_wave = sound_speed./freq; %wavelength
k = 2*pi./(lambda_wave); %wavenumber
t = 0; %unnecessary phase shift parameter

generalAngle = true; %a 'general angle' in this regard refers to that the
% direction of the sound propagation is, in the algorithm, assumed to be 
% the same everywhere, while the simulated sound is defined as coming from
% a point sound source and is thus incorrect; have yet to properly
% implement any improvement to this though, keep it as 'true'
%%

if db2pow(pressureNoise) == 0
    noise_rand = false;
else
    noise_rand = true;
end


for freqIdx = 1:size(freq,2)

    for iterIdx = 1:nIter
        %position of the (point) sound source
        srcpos = randomPointSourcePos(10,2,nSrc);
        %the 'measured' position of the sound source, with an error in msrmnt
        srcpos_noise = sourcePosNoise(srcpos, noiseError);

        %similar to above, mic coords and measured mic coords w/ error
        [mic_coords_real, mic_idx] = findRandMicsROI(allcoords, nmic, ROI_coords);
        mic_coords_measured = mic_coords_real + normrnd(0, mic_stdev ...
            , size(mic_coords_real));


        %calculating the K-matrix
        [KDir, KNoDir, KGauss] = kmatDir(mic_coords_measured, k(freqIdx), srcpos_noise, beta);

        %calculating the kappa vector
        [kappaDir, kappaNonDir, kappaGauss] = kernelFunc(...
            allcoords, mic_coords_measured, srcpos_noise, k(freqIdx), beta, generalAngle);

        %The real sound field, P, which is unknown
        P = simRealSoundField(allcoords,srcpos,A,k(freqIdx),w(freqIdx),t);

        %The pressure recorded at the microphones
        P_mic = zeros(1,nmic);
        for srcIter = 1:size(srcpos,2)
            for micIter = 1:size(mic_coords_real,2)
                micDist = sqrt(sum((mic_coords_real(:,micIter)-srcpos(:,srcIter)).^2));

                %calculating the sound pressure at the exact position of
                %the microphones - not necessarily on the grid
                P_mic(micIter) = P_mic(micIter) + updateSField(A(srcIter), micDist, k(freqIdx), w(freqIdx), t);
            end
        end

        %if we chose to add noise, this adds noise to the recorded sound
        %pressure of the microphones. Adds noise relative to the rest of
        %the sound field, to achieve the desired SNR
        if noise_rand == true
            noise_var = mean(P.^2,'all')/db2pow(pressureNoise);
            for micIter2 = 1:size(mic_coords_real,2)
                P_mic(micIter2) = P_mic(micIter2) + normrnd(0,sqrt(noise_var));
            end
        end
        
        %interpolated pressure
        P_hhz    =   zeros(size(P));
        P_hhzDir =   zeros(size(P));
        P_Gauss  =   zeros(size(P));

        for srcIter = 1:size(srcpos,2)
            [pHhzDir, pHhz, pGauss] = interpAllGrid(...
                P_mic, nmic, lambda, ...
                KDir(:,:,srcIter), KNoDir, KGauss, kappaDir(:,:,:,:,srcIter), kappaNonDir, kappaGauss);
            P_hhz = P_hhz + pHhz;
            P_hhzDir = P_hhzDir + pHhzDir;
            P_Gauss = P_Gauss + pGauss;
        end

        % sum of all pressure values of all
        % validation microphones, ^2
        intReal = sum(abs((P.*ROI_coords).^2),'all');

        % sum of the difference between the real
        % pressure values and the interpolated values
        intErrorHhzDir      = sum(abs(((P-P_hhzDir).*ROI_coords).^2),'all');
        intErrorHhzNonDir   = sum(abs(((P-P_hhz).*ROI_coords).^2),'all');
        intErrorGauss       = sum(abs(((P-P_Gauss).*ROI_coords).^2),'all');

        % SDR between the two
        SDR_hhzDir         = db(intReal/intErrorHhzDir, 'power');
        SDR_hhzNonDir      = db(intReal/intErrorHhzNonDir, 'power');
        SDR_Gauss          = db(intReal/intErrorGauss, 'power');


        SDR(freqIdx, iterIdx, :) = [SDR_hhzDir, SDR_hhzNonDir, SDR_Gauss];

    end
end
%%


figure(17)
hold off
plot(freq,squeeze(mean(SDR,2)), '--o')
hold on
plot([freq(1), freq(end)], [0, 0], '--')
beta_plot = beta;
xlabel('Frequency')
ylabel('SINAD (dB)')
Legends = compose('Directional Helmholtz kernel, \beta = %d', beta_plot);
legend({'Directional Helmholtz (\beta = 6)', 'Helmholtz','Gaussian', 'No interpolation'})
hold off

%%

z_plot = 25;
x = allcoords(:,:,:,1);
y = allcoords(:,:,:,2) ; 
plot_field(5004, 'Real sound field, at the mean z value of ROI', x(:,:,z_plot), y(:,:,z_plot), real(P(:,:,z_plot)), 0);
plot_field(5005, 'Real sound field, at the mean z value of ROI, interp', x(:,:,z_plot), y(:,:,z_plot), real(P_hhzDir(:,:,z_plot)), 0);

