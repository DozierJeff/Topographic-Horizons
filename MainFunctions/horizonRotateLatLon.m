function [SForward,SBackward] = horizonRotateLatLon(angleToRotate,Z,R,useParallel,varargin)
% [SForward,SBackward] = horizonRotateLatLon(angleToRotate,Z,R,useParallel [,E])
%
% Horizons in direction after grid rotatated by angleToRotate, for elevation
% grids in geographic format, as specified by R, with uniform angular spacing
% in latitude and uniform angular spacing in longitude (but spacing for
% latitude and longitude can be different).
% Values returned are the horizon angles and the azimuths for each point
% (which can differ from some analytic solution based on rotation angle).
%
% To pick a rotation angle based on the desired azimuth, use the
% rotationAngleFromAzimuth function.
%
% Input
%   angleToRotate - angle to rotate the grid
%   Z - elevation, either single or double
%   R - raster reference
%   useParallel - either true or false or 1 or 0, parallelizes the analyses
%       along the columns of the rotated elevation grid
% Optional input, Name-Value pairs
%   'E' - reference ellipsoid (output from referenceEllipsoid), default is
%       'WGS84', otherwise needed only if MATLAB version is earlier than R2020b
%       because from R2020b onward, this information is in the R.GeographicCRS
%       field
%   'method' - interpolation method for rotating the grid, 'nearest'
%       (default) or 'bilinear'
%
% Output, forward and backward structures with
%   sin of horizon angles
%   mean azimuth of angles
%   distances to horizons
%   ("forward" is defined as the angle above which the points are
%   hidden when facing the sun at its azimuth)

p = inputParser;
addRequired(p,'angleToRotate',@(x) isnumeric(x) && isscalar(x) &&...
    abs(x)<=180)
addRequired(p,'Z',@(x) isnumeric(x) && ismatrix(x) && isfloat(x))
addRequired(p,'R',@(x) contains(class(x),'rasterref') &&...
    contains(class(R),'geographic','IgnoreCase',true))
addRequired(p,'useParallel',@(x) isscalar(x) &&...
    (islogical(x) || isnumeric(x)))
addParameter(p,'E',referenceEllipsoid('wgs84'),...
    @(x) contains(class(x),'referenceEllipsoid') ||...
    (isnumeric(x) && length(x)==2))
addParameter(p,'method','nearest',@ischar)
parse(p,angleToRotate,Z,R,useParallel,varargin{:})

assert(nnz(isnan(Z))<numel(Z),'oops, the whole Z matrix is NaN');

angleToRotate = p.Results.angleToRotate;
Z = p.Results.Z;
R = p.Results.R;
useParallel = logical(p.Results.useParallel);
method = p.Results.method;
assert(strcmpi(method,'nearest'),...
    '''method'' ''nearest'' is currently the only one supported, ''bilinear'' is problematic')
if verLessThan('map','5.0')
    E = p.Results.E;
    if isnumeric(E)
        refEllipsoid = oblateSpheroid;
        refEllipsoid.SemimajorAxis = E(1);
        refEllipsoid.Eccentricity = E(2);
    else
        refEllipsoid = E;
    end
else
    refEllipsoid = R.GeographicCRS.Spheroid;
end

% latitudes and longitudes of the input grid
if verLessThan('map','5.1')
    [xIntrinsic,yIntrinsic] = meshgrid(1:size(Z,2),1:size(Z,1));
    [lat,lon] = intrinsicToGeographic(R,xIntrinsic,yIntrinsic);
else
    [lat,lon] = geographicGrid(R);
end

if mod(angleToRotate,90)==0
    [SForward,SBackward] = internal90(angleToRotate,lat,lon,Z,...
        refEllipsoid,useParallel);
else
    [tmpS] = internalHorizonPts(angleToRotate,lat,lon,Z,...
        refEllipsoid,useParallel,method);
    % cleanup rotated back to original grid
    [SForward,SBackward] = cleanupRotated(angleToRotate,tmpS,R,method);
end
end