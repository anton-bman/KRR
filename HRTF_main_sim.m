clear
clc



showPlot = false;
randomHRTF = false;
%the srcpos measured here is not used at all
[micCoords, earpos_left, earpos_right, ~] = getRealCoordinates2(showPlot);
valCoords = [earpos_left earpos_right];

%micCoords = [micCoords -micCoords]
B = [0; 0.05; 0];
A = earpos_left;

%we find the rotation vector that moves our origin-centered positions of
%the ears to [0 dist/2 0], where dist is the distance between the ears, and
%we use it to rotate all the coordinates.
rotMat = vrrotvec2mat(vrrotvec(A,B));

newValCoords = rotMat*valCoords;
newMicCoords = rotMat*micCoords;

newEarPosLeft = rotMat*earpos_left;
newEarPosRight = rotMat*earpos_right;

micCoords = newMicCoords;
%valCoords = newValCoords;
valCoords = [0 0; 0.09 -0.09; 0 0];

origin_point = [0 0 0]';



nFreqs = 1; %number of frequencies used in the signal
recLength = 5; %time of the recording used


%%%~~~~~~~------ frequencies used etc -----~~~~~~~%%%
sound_speed = 343; %roughly 343 at 20 degree C
freqs_rec = [100, 150, 200, 250, 300, 350, 400, 450, 500];%, 700, 900, 1100];

%%%~~~~~~~------ variables to define ------~~~~~~~%%%
beta = 3; %regularization parameter for the direction
lambda = 10^(-2); %regularization parameter
generalAngle = true;

%%%~~~~~~~----- Angle calculation etc -----~~~~~~~%%%

nIter = 20;

%%%~~~~~~~------ SOFA/HRTF stuff  etc -----~~~~~~~%%%
SOFAstart;
% File name of the RIEC file
RIECfn  = 'RIEC_hrir_subject_080';
SOFAfile=[RIECfn '.sofa'];

% Loading the full object
disp(['Loading full object: ' SOFAfile]);
Obj=SOFAload(SOFAfile);


% File name of the RIEC file

RIECfn2  = 'RIEC_hrir_subject_080';
SOFAfile2=[RIECfn '.sofa'];
disp(['Loading full object: ' SOFAfile2]);
Obj2=SOFAload(SOFAfile2);


%Obj.SourcePosition(:,1) = round(Obj.SourcePosition(:,1));
%Obj.SourcePosition(:,2) = round(Obj.SourcePosition(:,2));

%dist1 = mean(Obj.SourcePosition(:,3));

Fs = Obj.Data.SamplingRate; %samples per second
delayHRTF = 108; %delay in the HRTF, post-processing etc....

for iterIdxTot = 1:nIter
    if randomHRTF == true
        % File name of the RIEC file
        RIECfn2 = 'RIEC_hrir_subject_0';
        subject_id = randi(20,1);
        if subject_id < 10
            RIECfn2 = [RIECfn2 '0'];
        end
        SOFAfile2=[RIECfn2 int2str(subject_id) '.sofa'];
        disp(['Loading full object: ' SOFAfile2]);
        Obj2=SOFAload(SOFAfile2);
    end
    idx = iterIdxTot%randi(size(Obj2.SourcePosition(:,2),1),1);
    [srcpos(1,1),srcpos(2,1),srcpos(3,1)] = sph2cart(Obj.SourcePosition(idx,1),Obj.SourcePosition(idx,2),Obj.SourcePosition(idx,3));
    %srcpos
    dist1 = Obj.SourcePosition(idx,3);
    dist2 = Obj2.SourcePosition(idx,3);
    for freqIdxTot = 1:size(freqs_rec,2)
        
        %%%~~~~~~~------ SOFA/HRTF stuff  etc -----~~~~~~~%%%
        % Get index of measurements with the same directions
        %     idx=find(Obj.SourcePosition(:,1)==roundAzi(1) & Obj.SourcePosition(:,2)==roundEle(1));
        %    idx2=find(Obj.SourcePosition(:,1)==roundAzi(2) & Obj.SourcePosition(:,2)==roundEle(2));
        %     [~,idx] = min( sqrt( (Obj.SourcePosition(:,1)-roundAzi).^2 + (Obj.SourcePosition(:,2)-roundEle).^2 ) )
        %     [~,idx2] = min( sqrt( (Obj2.SourcePosition(:,1)-roundAzi).^2 + (Obj2.SourcePosition(:,2)-roundEle).^2 ) )

        % Extract and plot the fully loaded data
        RIEC_hrirL=squeeze(Obj.Data.IR(idx,1,:)); %left ear
        RIEC_hrirR=squeeze(Obj.Data.IR(idx,2,:)); %right ear

        % Extract and plot the fully loaded data
        RIEC_hrirL_calc=squeeze(Obj2.Data.IR(idx,1,:)); %left ear
        RIEC_hrirR_calc=squeeze(Obj2.Data.IR(idx,2,:)); %right ear

        %%% sound source that is being sent out at srcpos
        %     [primary_sources,Fs] = audioread('281_383_661_freq_combined.wav');
        %freqs = randi(500,1,3)+100;
        %freqs = 900;
        freqs = freqs_rec(freqIdxTot);
        [primary_sources,Fs] = newFileFreq(freqs,Fs, recLength);


        %downsamplingRatio = 2;
        %primary_sources = downsample(primary_sources,downsamplingRatio);
        %Fs = Fs./downsamplingRatio;
        %%%
        clear micCoords_pressure_timedomain;
        clear micCoords_pressure_timedomain_val;
        clear micCoords_pressure_timedomain_val_noear;

        for micIdx2 = 1:size(micCoords,2)
            r = sqrt(sum((srcpos-micCoords(:,micIdx2)).^2, 'all'));
            delay = r/(sound_speed)*Fs;
            transferfunction = zeros(1,ceil(max(delay)));
            transferfunction(ceil(max(delay))) = 1/r;

            %the 'real' sound pressure in each point, assuming the body won't change
            %the sound in any way aside from in the ear
            micCoords_pressure_timedomain(:,micIdx2) = filter(transferfunction,1,primary_sources);
        end

        for micValIdx2 = 1:size(valCoords,2)
            r_valMicDist(micValIdx2) = sqrt(sum((srcpos).^2, 'all'));
            delay = r_valMicDist(micValIdx2)/(sound_speed)*Fs;
            transferfunction = zeros(1,ceil(max(delay)));
            transferfunction(ceil(max(delay))) = 1/r_valMicDist(micValIdx2);
            %micCoords_pressure_timedomain_val(:,micValIdx2) = filter(transferfunction,1,primary_sources);

            if micValIdx2 == 1
                RIEC_temp = RIEC_hrirL;
            else
                RIEC_temp = RIEC_hrirR;
            end

            
            %the real sound in the position of the ear if the ear/head was not there
            micCoords_pressure_timedomain_val_noear(:,micValIdx2) = filter(transferfunction,1,primary_sources);

            %the real sound in the ear (simulated w/ HRIR)
            temp = conv(micCoords_pressure_timedomain_val_noear(:,micValIdx2),RIEC_temp);
            micCoords_pressure_timedomain_val(:,micValIdx2) = temp((1+delayHRTF):(size(primary_sources,1)+delayHRTF));

        end



        %miccoord pressure, validation miccoord pressure, primary sound source
        %pressure

        %     %time, when recording "starts"
        %     startTime = 1;
        %     primary_sources = primary_sources(startTime*Fs+1:end,:);
        %     micCoords_pressure_timedomain = micCoords_pressure_timedomain(round(startTime*Fs)+1:end,:);
        %     micCoords_pressure_timedomain_val = micCoords_pressure_timedomain_val(round(startTime*Fs)+1:end,:);
        %     micCoords_pressure_timedomain_val_noear = micCoords_pressure_timedomain_val_noear(round(startTime*Fs)+1:end,:);


        %%

        %     %we know there's only 3 frequencies
        %     nFreqs = 3;


        %we assume that every recording has the same frequencies, or we need to
        %loop through it all
        P = size(primary_sources,1);
        %N-point fft, padded with zeros if prim_sources2 has less than
        % N points and truncated if it has more.
        %Y   = fftshift( abs(fft(prim_sources))/size(prim_sources,1) ).^2; %ampl
        ff  = (0:P-1)/P;                    % Frequency grid
        Y   = ( abs(fft(micCoords_pressure_timedomain(:,1),P))/P ).^2; %ampl
        Y2  = (    (fft(micCoords_pressure_timedomain(:,1),P)) ); %ampl


        %finding the frequencies
        for i = 1:size(primary_sources,2)
            %replace this with any function that finds the largest peaks in
            %a dataset
            temp_peak = findpeaks2(Y(:,i),nFreqs*2);
            fPER(:,i) = temp_peak((1:2:nFreqs*2)+1); fPER(:,i)  = ff(fPER(:,i));
            for iterIdx = 1:size(fPER,1)
                tempIdx(iterIdx) = find(ff == fPER(iterIdx,i));
                freq(i,iterIdx) = fPER(iterIdx,i)*Fs;
            end
        end

        for freqIdx = 1:size(freq,2)


            lambda_wave = sound_speed/freq(freqIdx); %wavelength
            k = 2*pi/(lambda_wave); %wavenumber
            w = 2*pi*freq(freqIdx); %angular frequency


            micCoords_noise = micCoords;% + normrnd(0, mic_stdev*sqrt(pi/6), size(micCoords))
            %calculating the K-matrix
            [KDir, KNoDir, KGauss] =  kmatDir2(micCoords_noise , k, srcpos, beta, generalAngle);

            [kappaDir, kappaNonDir, kappaGauss] = kernelFunc( ...
                origin_point, micCoords_noise, srcpos, k, beta, generalAngle);

            [kappaDirAll, kappa2, kappa3] = kernelFunc( ...
                valCoords, micCoords_noise, srcpos, k, beta, generalAngle);
            %     [kappaDir_indiv, kappaNonDir_indiv, kappaGauss_indiv] = kernelFunc( ...
            %         valMicCoords, micCoords, srcpos, k, beta, false);

            for micIdx2 = 1:size(micCoords,2)
                micCoords_pressure(micIdx2) = goertzel(micCoords_pressure_timedomain(:,micIdx2), round(freq(freqIdx)/Fs*P)+1);
            end


            [pHhzDir(:,:,freqIdx,freqIdxTot), pHhz(:,:,freqIdx,freqIdxTot), pGauss(:,:,freqIdx,freqIdxTot)] = interpAll(micCoords_pressure, ...
                origin_point, lambda, ...
                KDir, KNoDir, KGauss, kappaDir, ...
                kappaNonDir, kappaGauss);

            [interpToEarCoords(:,:,freqIdx,freqIdxTot), ~, ~] = interpAll(micCoords_pressure, ...
                valCoords, lambda, ...
                KDir, KNoDir, KGauss, kappaDirAll, ...
                kappa2, kappa3);

        end

        %%

        interpSoundHhzDir = zeros(size(primary_sources,1), size(valCoords,2));
        interpSoundHhz = zeros(size(primary_sources,1), size(valCoords,2));
        interpSoundGauss = zeros(size(primary_sources,1), size(valCoords,2));

        testInterpSoundHhzDir = zeros(size(primary_sources,1), size(valCoords,2));

        for freqIdx = 1:nFreqs
            tempPressureHhzDir = pHhzDir(:,:,freqIdx, freqIdxTot);
            tempPressureHhz = pHhz(:,:,freqIdx, freqIdxTot);
            tempPressureGauss = pGauss(:,:,freqIdx, freqIdxTot);

            tempPressureEarPos = interpToEarCoords(:,:,freqIdx,freqIdxTot);
            for valMicIdx = 1:size(valCoords,2)

                %             temp = fft(micCoords_pressure_timedomain_val_noear(:,valMicIdx));
                %directional Hhz kernel
                freqPowerHhzDir = zeros(size(primary_sources,1),1);
                freqPowerHhzDir(round(freq(freqIdx)/Fs*P)+1) = tempPressureHhzDir;
                freqPowerHhzDir(size(freqPowerHhzDir,1) - round(freq(freqIdx)/Fs*P) + 1) = tempPressureHhzDir';

                %Hhz kernel
                freqPowerHhz = zeros(size(primary_sources,1),1);
                freqPowerHhz(round(freq(freqIdx)/Fs*P)+1) = tempPressureHhz;
                freqPowerHhz(size(freqPowerHhz,1) - round(freq(freqIdx)/Fs*P) + 1) = tempPressureHhz';

                %Gauss kernel
                freqPowerGauss = zeros(size(primary_sources,1),1);
                freqPowerGauss(round(freq(freqIdx)/Fs*P)+1) = tempPressureGauss;
                freqPowerGauss(size(freqPowerGauss,1) - round(freq(freqIdx)/Fs*P) + 1) = tempPressureGauss';

                %DirHhz kernel without the HRTF stuff
                testPower_dirHhz = zeros(size(primary_sources,1),1);
                testPower_dirHhz(round(freq(freqIdx)/Fs*P)+1) = tempPressureEarPos(valMicIdx);
                testPower_dirHhz(size(freqPowerHhzDir,1) - round(freq(freqIdx)/Fs*P) + 1) = tempPressureEarPos(valMicIdx)';


                if valMicIdx == 1
                    RIEC_temp = RIEC_hrirL;
                else
                    RIEC_temp = RIEC_hrirR;
                end

                soundPressureTimeDirHhz = ifft(freqPowerHhzDir);
                soundPressureTimeHhz = ifft(freqPowerHhz);
                soundPressureTimeGauss = ifft(freqPowerGauss);

                temp = conv(soundPressureTimeDirHhz,RIEC_temp);
                soundPressureTimeDirHhz = temp((1+delayHRTF):(size(primary_sources,1)+delayHRTF));

                temp = conv(soundPressureTimeHhz,RIEC_temp);
                soundPressureTimeHhz = temp((1+delayHRTF):(size(primary_sources,1)+delayHRTF));

                temp = conv(soundPressureTimeGauss,RIEC_temp);
                soundPressureTimeGauss = temp((1+delayHRTF):(size(primary_sources,1)+delayHRTF));




                testPressureTimeDirHhz = ifft(testPower_dirHhz);
                %final interpolation
                interpSoundHhzDir(:,valMicIdx) = interpSoundHhzDir(:,valMicIdx) + soundPressureTimeDirHhz;
                interpSoundHhz(:,valMicIdx) = interpSoundHhz(:,valMicIdx) + soundPressureTimeHhz;
                interpSoundGauss(:,valMicIdx) = interpSoundGauss(:,valMicIdx) + soundPressureTimeGauss;

                testInterpSoundHhzDir(:,valMicIdx) = testInterpSoundHhzDir(:,valMicIdx) + testPressureTimeDirHhz;
            end



            
        end

        %remove first 1000 samples to get rid of convolution and filtering
        %delay stuff


        testInterpSoundHhzDir = testInterpSoundHhzDir((1000+1):(end-400),:);
        interpSoundHhzDir = interpSoundHhzDir((1000+1):(end-400),:);
        interpSoundHhz = interpSoundHhz((1000+1:(end-400)),:);
        interpSoundGauss = interpSoundGauss((1000+1):(end-400),:);

        micCoords_pressure_timedomain_val = micCoords_pressure_timedomain_val((1000+1):(end-400),:);
        %     figure(3338)
        %     hold off
        %     plot(micCoords_pressure_timedomain_val(:,1))
        %     hold on
        %     plot(micCoords_pressure_timedomain_val(:,1) - interpSoundGauss(:,1))
        %     plot(micCoords_pressure_timedomain_val(:,1) - interpSoundHhz(:,1))
        %     plot(micCoords_pressure_timedomain_val(:,1) - interpSoundHhzDir(:,1))

        %%
        intReal = sum((micCoords_pressure_timedomain_val).^2,1);

        intErrorHhzDir = sum((micCoords_pressure_timedomain_val - interpSoundHhzDir).^2,1);
        intErrorHhz = sum((micCoords_pressure_timedomain_val - interpSoundHhz).^2,1);
        intErrorGauss = sum((micCoords_pressure_timedomain_val - interpSoundGauss).^2,1);

        % SNR for free field
        intErrorHhzDir_freeField = sum((micCoords_pressure_timedomain_val - testInterpSoundHhzDir).^2,1);


        SNR_hhzDir(freqIdxTot,iterIdxTot,:) = db(intReal./intErrorHhzDir, 'power');
        SNR_hhz(freqIdxTot,iterIdxTot,:) = db(intReal./intErrorHhz, 'power');
        SNR_Gauss(freqIdxTot,iterIdxTot,:) = db(intReal./intErrorGauss, 'power');

        SNR_hhzDir_freeField(freqIdxTot,iterIdxTot,:) = db(intReal./intErrorHhzDir_freeField, 'power');


        figure(123123)
        hold off
        plot(micCoords_pressure_timedomain_val(:,1))
        hold on
        plot(micCoords_pressure_timedomain_val(:,1)-interpSoundGauss(:,1))
        plot(micCoords_pressure_timedomain_val(:,1)-interpSoundHhz(:,1))
        plot(micCoords_pressure_timedomain_val(:,1)-interpSoundHhzDir(:,1))
        hold off
        savedFreq(iterIdxTot,:) = freqs;
    end
end

%%

hhzDir_SNR = mean(mean(SNR_hhzDir))
hhz_SNR = mean(mean(SNR_hhz))
gauss_SNR = mean(mean(SNR_Gauss))


mean_SDR = [mean(mean(SNR_hhzDir,2),3) mean(mean(SNR_hhz,2),3) mean(mean(SNR_Gauss,2),3) mean(mean(SNR_hhzDir_freeField,2),3)];

%%
% figure(1338)
%plot(micCoords_pressure_timedomain_val_noear(:,1))
%plot(micCoords_pressure_timedomain_val(:,1))
% hold on
% plot(testInterpSoundHhzDir(:,1))
% plot(interpSoundHhzDir(:,1))
% hold off

%%
figure(15)
hold off
plot(freqs_rec(1:9),mean_SDR(1:9,:), '--o')
hold on
plot([freqs_rec(1),freqs_rec(9)], [0, 0], '--')
beta_plot = [beta, 0];
xlabel('Frequency')
ylabel('SINAD (dB)')
Legends = compose('Proposed method (beta = %d)', beta_plot);
legend({char(Legends(1)), char(Legends(2)),'Proposed method (Gauss)','KRR, without HRTF (beta = 3)', 'No interpolation'})
hold off
% figure(13)
% hold on
% plot(freqs_rec,mean_SDR(:,1))
% beta_plot = [15, 0];
% hold off
% Legends = compose('Helmholtz kernel, beta = %d', beta_plot);
% newlegend = compose('Helmholtz kernel, beta = %d', beta);
% legend({char(Legends(1)), char(Legends(2)),'Gaussian Kernel','Free field estimate', 'no interpolation', char(newlegend(1)) })

%%

legend({'Proposed method (\beta = 3, HRTF)', 'Proposed method (\beta = 0, HRTF)','Proposed method (Gauss, HRTF)','KRR (\beta = 3, no HRTF)', 'No interpolation'},'Location','northeast');
