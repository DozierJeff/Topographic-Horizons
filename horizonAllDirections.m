function [A,H,varargout] = horizonAllDirections(Z,R,varargin)
% [azm,horzAng] = horizonAllDirections(Z,R [,name-value pairs])
% [azm,horzAng,horzDis] = horizonAllDirections(Z,R [,name-value pairs])
%calculates horizons for elevation grid in all azimuths around circle,
%either -180 to +180 degrees counter-clockwise or 0 to 360 degrees clockwise
%
%Values returned in the structure SH are the horizon angles and the azimuths
%for each point, along with the distances (meters) from each point to its
%horizon.
%
%Input
%   Z - elevation grid in meters, single or double precision
%   R - raster reference for the elevation grid, can be a MapCellsReference
%       (or MapPostingsReference) object, in which case the grid is projected
%       and distances between points are based on the projected coordinates,
%       or a GeographicCellsReference (or GeographicPostingsReference object),
%       in which case the distances between points are calculated from the
%       great-circle distances
%       If on MATLAB version R2020b or later, the R object should include
%       the ProjectedCRS or GeographicCRS field but if omitted or empty the
%       'proj' or 'planet' optional input can be used
%Optional input, name-value pairs, case-insensitive
%   'proj' - projection structure or projcrs object, needed if and only if
%       R is a MapCellsReference or MapPostingsReference object without the
%       non-empty ProjectedCRS field
%   'planet' - applicable if and only if R is a GeographicCellsReference or
%       GeographicPostingsReference object and does not contain a non-empty
%       GeographicCRS field, default is 'earth' with a WGS84 ellipsoid
%   'nHorz' - number of horizon directions over full circle, default 64,
%       but adjusted to multiple of 4 otherwise
%   'selectAng' - instead of full circle, vector of rotation angles to
%       evaluate, overrides 'nHorz' if provided
%   'parallel' - followed by 'profile' to parallelize along the columns of
%       the rotated grid or by 'rotate' to run the suite of rotation angles
%       in parallel, default '' in which case no parallel processing is
%       employed
%   'method' - argument for interpolation in the rotation of the grid,
%       default 'nearest', other option in a future implementation might be
%       'bilinear', but unintended effects cause it not to be implemented now
%   'verbose' - if true, indicates when each rotation is done, default
%       false
% Output
%   azm - vector of azimuths, 1 for each rotation angle in nHorz directions
%   horzAng - 3-D horizons in BSQ (band-sequential) format, the 3rd dimension
%       specifying the elevation angle to the horizon (i.e. elevation angle
%       above the horizontal)
% Optional output
%   horzDis - distances to the horizons in meters corresponding to output H
%
% NOTE
% The 3D outputs are in 'bsq' (band-sequential) format, 2D locations stacked
% on horizons for the set of azimuths. You can use bsq2bip and bip2bsq to
% switch between band-sequential and band-interleaved-by-pixel formats.
% The way azimuths are handled assumes that the origin of the grid is the
% NW corner. Depending on the return value from azimuthPreference, the
% azimuths either range over ±180° counter-clockwise with 0° south, or over
% 0° to 360° clockwise with 0° north.

% flag to return distances
returnD = nargout>2;

% defaults
defaultH = 64;
defaultPlanet = 'wgs84';
defaultInterp = 'nearest';
minargs = 2;
maxargs = 16;
narginchk(minargs,maxargs)

p = inputParser;
addRequired(p,'Z',@(x) isnumeric(x) && ismatrix(x) && isfloat(x))
addRequired(p,'R',@(x) contains(class(x),'rasterref'))
addParameter(p,'nhorz',defaultH,@(x) isnumeric(x) && isscalar(x) && x>0)
addParameter(p,'planet',defaultPlanet,@ischar)
addParameter(p,'proj',struct([]),@(x) isstruct(x) ||...
    contains(class(x),'projcrs'))
addParameter(p,'parallel','',@ischar)
addParameter(p,'method',defaultInterp,@ischar)
addParameter(p,'verbose',false,@(x) islogical(x) || isnumeric(x))
addParameter(p,'selectang',[],@(x) isnumeric(x) && isvector(x))
parse(p,Z,R,varargin{:});

E = referenceEllipsoid(p.Results.planet);

% set parallel option
useParallelAlongProfiles = false;
useParallelAcrossRotations = false;
if strncmpi(p.Results.parallel,'pro',3)
    useParallelAlongProfiles = true;
elseif strncmpi(p.Results.parallel,'rot',3)
    useParallelAcrossRotations = true;
elseif ~isempty(p.Results.parallel)
    error('''parallel'' option ''%s'' not recognized',p.Results.parallel)
end

% verbose
verbose = logical(p.Results.verbose);

%other inputs
nHorz = p.Results.nhorz;
method = p.Results.method;
assert(strcmpi(method,'nearest'),...
    '''method'' ''nearest'' is currently the only one supported, ''bilinear'' is problematic')
if contains(class(R),'geographic','IgnoreCase',true)
    useLatLon = true;
    if ~isempty(p.Results.proj)
        warning('R (2nd argument) is of class %s, so ''proj'' should be empty, it''s ignored',class(R))
    end
    proj = struct([]);
else
    useLatLon = false;
    assert(contains(class(R),'map.rasterref.Map','IgnoreCase',true),...
        'R (2nd argument) must be a Map or Geographic raster reference')
    if verLessThan('map','5.0')
        assert(~isempty(p.Results.proj),...
            '''proj'' can''t be empty for MATLAB versions before R2020b')
        proj = p.Results.proj;
    else
        proj = R.ProjectedCRS;
    end
end

% nHorz must be multiple of 4
if isempty(p.Results.selectang)
    if nHorz<4
        nHorz = 4;
        warning('''nHorz'' must be >=4, set to 4')
    elseif mod(nHorz,4)~=0
        nHorz = 4*ceil(nHorz/4);
        warning('''nHorz'' must be a multiple of 4, set to %d',nHorz)
    end

    % rotation angles, need just southern half of directions because horz2d
    % computes forward and backward
    ang = linspace(-90,90,1+nHorz/2);
    ang(1) = []; % 90 and -90 gotten on same rotation
else
    ang = p.Results.selectang;
    nHorz = length(ang);
end

% 3D arays to hold horizons in forward and backward directions for each azimuth
HgridF = zeros(size(Z,1),size(Z,2),length(ang),'single');
HgridB = zeros(size(HgridF),'like',HgridF);
AvecF = zeros(length(nHorz),1);
AvecB = zeros(length(nHorz),1);
if returnD
    DgridF = zeros(size(HgridF),'like',HgridF);
    DgridB = zeros(size(HgridF),'like',HgridF);
end

% compute, depending on parallel option
if useParallelAcrossRotations
    % useParallelAlongProfiles passed to the rotational parfor loops, even
    % though this variable is set to false when useParallelAcrossRotations
    % is true, but could be true in situations where nested parallelism
    % could be implemented
    if useLatLon
        if verLessThan('map','5.0')
            passVals = {Z,R,useParallelAlongProfiles,method,verbose,E};
        else
            passVals = {Z,R,useParallelAlongProfiles,method,verbose};
        end
    else
        if verLessThan('map','5.0')
            passVals = {Z,R,useParallelAlongProfiles,method,verbose,proj};
        else
            passVals = {Z,R,useParallelAlongProfiles,method,verbose};
        end
    end
    W = parallel.pool.Constant(passVals);
    if useLatLon
        parfor h=1:length(ang)
            pv = W.Value;
            if verLessThan('map','5.0')
                [SF,SB] = horizonRotateLatLon(ang(h),pv{1},pv{2},pv{3},...
                    'method',pv{4},'E',pv{6},'verbose',pv{5});
            else
                [SF,SB] = horizonRotateLatLon(ang(h),pv{1},pv{2},pv{3},...
                    'method',pv{4},'verbose',pv{5});
            end

            HgridF(:,:,h) = single(SF.horzAng);
            HgridB(:,:,h) = single(SB.horzAng);
            if returnD
                DgridF(:,:,h) = single(SF.horzDis);
                DgridB(:,:,h) = single(SB.horzDis);
            end
            AvecF(h) = SF.azm;
            AvecB(h) = SB.azm;
        end
    else
        parfor h=1:length(ang)
            pv = W.Value;
            if verLessThan('map','5.0')
                [SF,SB] = horizonRotateProj(ang(h),pv{1},pv{2},pv{3},...
                    'proj',pv{6},'method',pv{4},'verbose',pv{5});
            else
                [SF,SB] = horizonRotateProj(ang(h),pv{1},pv{2},pv{3},...
                    'method',pv{4},'verbose',pv{5});
            end
            HgridF(:,:,h) = single(SF.horzAng);
            HgridB(:,:,h) = single(SB.horzAng);
            if returnD
                DgridF(:,:,h) = single(SF.horzDis);
                DgridB(:,:,h) = single(SB.horzDis);
            end
            AvecF(h) = SF.azm;
            AvecB(h) = SB.azm;
        end
    end
else
    if useLatLon
        for h=1:length(ang)
            if verLessThan('map','5.0')
                [SF,SB] = horizonRotateLatLon(ang(h),Z,R,useParallelAlongProfiles,...
                    'E',E,'method',method,'verbose',verbose);
            else
                [SF,SB] = horizonRotateLatLon(ang(h),Z,R,useParallelAlongProfiles,...
                    'method',method,'verbose',verbose);
            end
            HgridF(:,:,h) = single(SF.horzAng);
            HgridB(:,:,h) = single(SB.horzAng);
            if returnD
                DgridF(:,:,h) = single(SF.horzDis);
                DgridB(:,:,h) = single(SB.horzDis);
            end
            AvecF(h) = SF.azm;
            AvecB(h) = SB.azm;
        end
    else
        for h=1:length(ang)
            if verLessThan('map','5.0')
                [SF,SB] = horizonRotateProj(ang(h),Z,R,useParallelAlongProfiles,...
                    'proj',proj,'method',method,'verbose',verbose);
            else
                [SF,SB] = horizonRotateProj(ang(h),Z,R,useParallelAlongProfiles,...
                    'method',method,'verbose',verbose);
            end
            HgridF(:,:,h) = single(SF.horzAng);
            HgridB(:,:,h) = single(SB.horzAng);
            if returnD
                DgridF(:,:,h) = single(SF.horzDis);
                DgridB(:,:,h) = single(SB.horzDis);
            end
            AvecF(h) = SF.azm;
            AvecB(h) = SB.azm;
        end
    end
end

% gather list of variables to clear
toClear = {'SF','SB','HgridF','HgridB'};
if returnD
    toClear = cat(2,toClear,{'DgridF','DgridB'});
end

%concatenate and sort by azimuths
Hcat = cat(3,HgridF,HgridB);
if returnD
    Dcat = cat(3,DgridF,DgridB);
end
Acat = cat(2,AvecF,AvecB);
[A,I] = sort(Acat);
H = zeros(size(Hcat),'like',Hcat);
toClear = cat(2,toClear,{'Hcat'});
if returnD
    toClear = cat(2,toClear,{'Dcat'});
    D = zeros(size(Dcat),'like',Dcat);
    for k=1:size(H,3)
        H(:,:,k) = Hcat(:,:,I(k));
        D(:,:,k) = Dcat(:,:,I(k));
    end
    mind = prctile(D(:),.5);
    D(D<mind) = mind;
else
    for k=1:size(H,3)
        H(:,:,k) = Hcat(:,:,I(k));
    end
end

% to prevent memory error, get rid of variables no longer needed
clear(toClear{:})

% smooth 3D H & D, making sure kernelSize is odd and is <= 7
kernelSize = round(max(3,min(7,length(A)/(128/5)))); % 3x3x3 minimum, 5x5x5 if nHorz=128
if mod(kernelSize,2)==0
    kernelSize = kernelSize+1;
end
H = smooth3(H,'gaussian',kernelSize);
if returnD
    D = smooth3(D,'gaussian',kernelSize);
end

% make periodic, place in structure, and return
[Anew,addFirst,addLast] = fillAzimuthEnds(A);
if length(Anew)==2
    pause
end
H = makeVariablePeriodic(addFirst,addLast,H,Anew);
if returnD
    D = makeVariablePeriodic(addFirst,addLast,D,Anew);
    varargout{1} = D;
end
A = Anew;

end

function [newA,addFirst,addLast] = fillAzimuthEnds(azm)
if min(azm)<0
    azrange = [-180 180];
else
    azrange = [0 360];
end
% add at beginning and/or end
addFirst = false;
addLast = false;
azm = azm(:).';
if min(azm)>min(azrange)
    azm = [min(azrange) azm];
    addFirst = true;
end
if max(azm)<max(azrange)
    azm = [azm max(azrange)];
    addLast = true;
end
newA = unique(azm); %should already, but just in case
end

function newX = makeVariablePeriodic(addFirst,addLast,X,azm)
% do nothing if already covers ends
if ~(addFirst || addLast)
    newX = X;
else
    X = bsq2bip(X);
    if length(azm)==size(X,1) % shouldn't, but just in case
        newX = bip2bsq(X);
    else
        newX = zeros(length(azm),size(X,2),size(X,3),'like',X);
        if addFirst && addLast
            newX(2:end-1,:,:) = X;
            newX(1,:,:) = (newX(2,:,:)+newX(end-1,:,:))/2;
            newX(end,:,:) = newX(1,:,:);
        elseif addFirst
            newX(2:end,:,:) = X;
            newX(1,:,:) = newX(end,:,:);
        else %addLast
            newX(1:end-1,:,:) = X;
            newX(end,:,:) = newX(1,:,:);
        end
        newX = bip2bsq(newX);
    end
end
end