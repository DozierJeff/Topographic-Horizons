function et = Demo_shade(Z,R,dt,useParallel)
%figure showing shading by slope and horizon for local datetime dt

tic; % start the timer

% location of sun
[declin,~,omega,~] = EarthEphemeris(dt);

% sun angle and azimuth at every point
if verLessThan('map','5.1')
[xI,yI] = meshgrid(1:R.RasterSize(2),1:R.RasterSize(1));
[lat,lon] = intrinsicToGeographic(R,xI,yI);
else
    [lat,lon] = geographicGrid(R);
end
[mu0,phi0] = sunang(lat,lon,declin,omega);
az = mean(phi0(:));
fprintf('This code %s reproduces Fig. 3\n',mfilename);
fprintf('solar zenith %f\n',acosd(mean(mu0(:))))
fprintf('solar azimuth %f counterclockwise from south\n',mean(phi0(:)))

% slope, aspect
[slope,aspect] = topographicSlope(double(Z),R);

% sun angles on slope
mu = sunslope(mu0,phi0,slope,aspect);

% horizons at azimuth
[angToRotate,~,~] = rotationAngleFromAzimuth(az,R);
[SForward,SBackward] = horizonRotateLatLon(angToRotate,double(Z),R,useParallel);
disp(SForward)
disp(SBackward)

% which in shade
shadeSlope = mu<=0;
shadeHorizon = sind(SForward.horzAng)>mu0;

% image, shaded relief, shade by slope, shade by horizon
% shaded relief
figure('Name','Fig. 3 Shaded by slope or horizon')
ax = setAxes(R,true); %#ok<NASGU>
X = zeros(size(Z));
X(shadeSlope) = 1;
X(shadeHorizon & ~shadeSlope) = 2;
Q = [242 242 242; 227 191 105; 77 51 204];
ax = setAxes(R,true); %#ok<NASGU>
geoshow(X,R,'DisplayType','Surface')
colormap(Q/255)
colorbar('Location','SouthOutside','Ticks',[0 1 2],'TickLabels',...
    {'no shade','slope','horizon'},'Box','off')

% percentages
N = numel(Z);
fprintf('fraction illuminated = %f\n',nnz(~(shadeSlope | shadeHorizon))/N);
fprintf('fraction shaded by slope = %f\n',nnz(shadeSlope)/N);
fprintf('fraction shaded by horizon and not by slope = %f\n',...
    nnz(shadeHorizon & ~shadeSlope)/N);

et = toc;

end