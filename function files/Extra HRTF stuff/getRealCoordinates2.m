function [mic_coords, earpos_left, earpos_right, srcpos] = getRealCoordinates2(showPlot)
% Gets the measured coordinates, and sets the origin as
% the mean coordinate of the in-ear microphones, instead of all microphones
% (as opposed to getRealCoordinates)
% showPlot = true will plot the room

[mic_coords, earpos_right,earpos_left,srcpos] = anton_calibrate_setup(showPlot);

offset = mean([earpos_right earpos_left],2);
mic_coords = mic_coords - offset;
earpos_right = earpos_right - offset;
earpos_left = earpos_left - offset;
srcpos = srcpos - offset;

end