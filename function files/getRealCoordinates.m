function [mic_coords, earpos_left, earpos_right, srcpos] = getRealCoordinates(showPlot)
% Gets the measured coordinates, and sets the origin as
% the mean coordinate of all microphones
% showPlot = true will plot the room

[mic_coords, earpos_right,earpos_left,srcpos] = anton_calibrate_setup(showPlot);

offset = mean(mic_coords,2);
mic_coords = mic_coords - offset;
earpos_right = earpos_right - offset;
earpos_left = earpos_left - offset;
srcpos = srcpos - offset;

end