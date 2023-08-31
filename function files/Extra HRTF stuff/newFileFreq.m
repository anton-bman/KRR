function [data,Fs] = newFileFreq(freqs, Fs, recLen)



amp = ones(1,size(freqs,2));
dt = 1/Fs;                   % seconds per sample
t = (0:dt:recLen-dt)';     % seconds
F1 = 100;                      % Sine wave frequency (hertz)
F2 = 250;
F3 = 550;

data = 0;
for freqIdx = 1:size(freqs,2)
    data = data + amp(freqIdx)*sin(2*pi*freqs(freqIdx)*t + 2*pi*rand);
end

data = data./max(data);