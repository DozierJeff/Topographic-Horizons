function S = internalHorizonPts(angToRotate,lat,lon,Z,refEllipsoid,UseParallel,method)
%rotation when not multiple of 90 degrees
assert(mod(angToRotate,90)~=0,...
    'rotation %f is a multiple of 90, use internal90 instead',angToRotate)
assert(strcmpi(method,'nearest'),'method must be ''nearest'', others not implemented')

%original NaN values in elevation
origNaN = isnan(Z);

%imrotate marks areas outsize the original grid as zeros, but data values
%could also be zeros. With nearest-neighbor rotation, we can set zero
%values in the data to a placeholder value (-Inf) and then reset after we
%have changed the values outside the original grid to NaN
lat(lat==0) = -Inf;
lon(lon==0) = -Inf;
Z(Z==0) = -Inf;
lat = imrotate(lat,angToRotate);
lon = imrotate(lon,angToRotate);
Z = imrotate(Z,angToRotate);
mask = lat==0;
lat(mask) = NaN;
lon(mask) = NaN;
Z(mask) = NaN;
lat(isinf(lat)) = 0;
lon(isinf(lon)) = 0;
Z(isinf(Z)) = 0;

% output arrays
holdH = nan(size(lat));
holdA = nan(size(lat,2),2); % for mean of azimuths along profile
holdDis = nan(size(lat));
backH = nan(size(holdH));
backA = nan(size(holdA));
backDis = nan(size(holdDis));

% minimum number of points in profile to enable calculation
minPts = 5;
if UseParallel
    parfor k=1:size(lat,2)
        [horzAng,horzDis,azm,npts] = horizonAlongProfile(lat(:,k),...
            lon(:,k),Z(:,k),refEllipsoid);
        if nnz(~isnan(horzAng))>=minPts
            holdH(:,k) = horzAng;
            holdA(k,:) = [atan2d(mean(sind(azm),'omitnan'),mean(cosd(azm),'omitnan'))...
                nnz(~isnan(npts))];
            holdDis(:,k) = horzDis;
            [horzAng,horzDis,azm,npts] = horizonAlongProfile(flipud(lat(:,k)),...
                flipud(lon(:,k)),flipud(Z(:,k)),refEllipsoid);
            backH(:,k) = horzAng;
            backA(k,:) = [atan2d(mean(sind(azm),'omitnan'),mean(cosd(azm),'omitnan'))...
                nnz(~isnan(npts))];
            backDis(:,k) = horzDis;
        end
    end
else
    for k=1:size(lat,2)
        [horzAng,horzDis,azm,npts] = horizonAlongProfile(lat(:,k),...
            lon(:,k),Z(:,k),refEllipsoid);
        if nnz(~isnan(horzAng))>=minPts
            holdH(:,k) = horzAng;
            holdA(k,:) = [atan2d(median(sind(azm),'omitnan'),...
                median(cosd(azm),'omitnan'))...
                nnz(~isnan(npts))];
            holdDis(:,k) = horzDis;
            [horzAng,horzDis,azm,npts] = horizonAlongProfile(flipud(lat(:,k)),...
                flipud(lon(:,k)),flipud(Z(:,k)),refEllipsoid);
            backH(:,k) = horzAng;
            backA(k,:) = [atan2d(median(sind(azm),'omitnan'),...
                median(cosd(azm),'omitnan'))...
                nnz(~isnan(npts))];
            backDis(:,k) = horzDis;
        end
    end
end
backH = flipud(backH);
backDis = flipud(backDis);
backA = flipud(backA);

% variables to return
S.origNaN = origNaN;
S.lat = lat;
S.lon = lon;
S.Forward.horzAng = holdH;
S.Forward.horzDis = holdDis;
S.Backward.horzAng = backH;
S.Backward.horzDis = backDis;
% weighted mean of the columwise azimuths
for p=1:2
    switch p
        case 1
            A = holdA;
        case 2
            A = backA;
    end
    if verLessThan('matlab','9.9')
        wt = A(:,2)./nansum(A,2); %#ok<*NANSUM>
        azmMean = atan2d(nansum(sind(A(:,1)).*wt),nansum(cosd(A(:,1)).*wt));
    else
        wt = A(:,2)./sum(A,2,'omitnan');
        azmMean = atan2d(sum(sind(A(:,1)).*wt,'omitnan'),...
            sum(cosd(A(:,1)).*wt,'omitnan'));
    end
    % atan2d returns +/- 180, fix if we use 0 to 360
    if ~azimuthPreference && azmMean<0
        azmMean = 360+azmMean;
    end
    switch p
        case 1
            S.Forward.azm = azmMean;
        case 2
            S.Backward.azm = azmMean;
    end
end
end