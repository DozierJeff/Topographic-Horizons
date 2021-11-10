function et = Demo_imageDEM(Z,R)
%images of DEM similar to Fig. 1

% intrinsic coordinates
if verLessThan('map','5.1')
    [xI,yI] = meshgrid(1:R.RasterSize(2),1:R.RasterSize(1));
    [lat,lon] = intrinsicToGeographic(R,xI,yI);
else
    [lat,lon] = geographicGrid(R);
end

tic; % start the timer

% image of the elevations
figure('Name','Fig. 1 left Elevation')
ax = setAxes(R,true); %#ok<NASGU>
geoshow(double(Z),R,'DisplayType','surface');
colormap(demcmap(double(Z)))
colorbar('SouthOutside')

% shaded relief
figure('Name','Fig.1 right Shaded Relief')
ax = setAxes(R,true); %#ok<NASGU>
surflsrm(lat,lon,double(Z),[45 315])

et = toc;

fprintf('This code %s reproduces the 2 images in Fig. 1\n',...
    mfilename)

end