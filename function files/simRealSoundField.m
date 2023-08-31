function P = simRealSoundField(allcoords,srcpos,A,k,w,t)
%simulates the real sound field, with source emitted from all the 
% positions in srcpos, with wave properties [A,k,w,t] corresponding to
% the amplitude, wave number, wave length and time, respectively.

%Pressure in the points of the sound field
Nall = size(allcoords);
Nx = Nall(1);
Ny = Nall(2);
Nz = Nall(3);

%pressure in the sound field
P = zeros(Nx,Ny,Nz);

xpos = srcpos(1,:);
ypos = srcpos(2,:);
zpos = srcpos(3,:);

x =  allcoords(:,:,:,1);
y =  allcoords(:,:,:,2);
z =  allcoords(:,:,:,3);


for j = 1:size(srcpos,2)
    %updating every point in the field, given the loudspeaker coordinates
    r = (((x-xpos(j)).^2) + (y-ypos(j)).^2 + (z-zpos(j)).^2).^(1/2);
    %P = P + update_sf3d(A(j),allcoords,xpos(j),ypos(j),zpos(j),k,w,t);

    
    P = P + updateSField(A(j),r, k, w, t);
end

end