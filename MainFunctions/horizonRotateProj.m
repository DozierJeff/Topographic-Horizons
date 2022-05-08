function [SForward,SBackward] = horizonRotateProj(angToRotate,Z,R,useParallel,varargin)
% [SForward,SBackward] = horizonRotateProj(angToRotate,Z,R,useParallel [,proj])
%
% Horizons in direction after grid rotatated by angToRotate, for elevation
% grids in a projection.
% Values returned are the horizon angles and the azimuths for each point
% (which can differ from some analytic solution based on rotation angle
% depending on projection).
%
% To pick a rotation angle based on the desired azimuth, use the
% rotationAngleFromAzimuth function
%
% Input
%   angToRotate - angle to rotate the grid
%   Z - elevation, either single or double
%   R - raster reference
%   useParallel - either true or false or 1 or 0, parallelizes the analyses
%       along the columns of the rotated elevation grid, default false
% Optional input, Name-Value pairs
%   'proj' - projection structure or projcrs object, needed only if MATLAB
%       version is earlier than R2020b because from R2020b onward, this
%       information is in the R.ProjectedCRS field
%   'method' - interpolation method for rotating the grid, 'nearest'
%       (default) or 'bilinear'
%   'verbose' - displays angle and azimuths of rotation during calculation
%       (default false)
%
% Output
%   Structures for the forward and backward directions, with fields
%   including:
%   elevation angles in degrees of forward and backward horizons
%   distances to the horizon for each point
%   azimuths to horizons, averaged over the grid

p = inputParser;
addRequired(p,'angToRotate',@(x) isnumeric(x) && isscalar(x) &&...
    abs(x)<=180)
addRequired(p,'Z',@(x) isnumeric(x) && ismatrix(x) && isfloat(x))
addRequired(p,'R',@(x) contains(class(x),'rasterref') &&...
    contains(class(x),'Map','IgnoreCase',true))
addRequired(p,'useParallel',@(x) isscalar(x) &&...
    (islogical(x) || isnumeric(x)))
addParameter(p,'proj',struct([]),@(x) contains(class(x),'projcrs') || isstruct(x))
addParameter(p,'method','nearest',@ischar)
addParameter(p,'verbose',false,@(x) islogical(x) || isnumeric(x))
parse(p,angToRotate,Z,R,useParallel,varargin{:})

assert(nnz(isnan(Z))<numel(Z),'oops, the whole Z matrix is NaN');
assert(strcmpi(p.Results.method,'nearest'),...
    '''method'' ''nearest'' is currently the only one supported, ''bilinear'' is problematic')

angToRotate = p.Results.angToRotate;
Z = p.Results.Z;
R = p.Results.R;
useParallel = logical(p.Results.useParallel);
verbose = logical(p.Results.verbose);
method = p.Results.method;

% latitudes and longitudes of the input grid
if verLessThan('map','5.1')
    [xIntrinsic,yIntrinsic] = meshgrid(1:size(Z,2),1:size(Z,1));
    [x,y] = intrinsicToWorld(R,xIntrinsic,yIntrinsic);
else
    [x,y] = worldGrid(R);
end
if verLessThan('map','5.0')
    proj = p.Results.proj;
    [lat,lon] = projinv(proj,x,y);
    E = proj.geoid;
    refEllipsoid = oblateSpheroid;
    refEllipsoid.SemimajorAxis = E(1);
    refEllipsoid.Eccentricity = E(2);
else
    [lat,lon] = projinv(R.ProjectedCRS,x,y);
    refEllipsoid = R.ProjectedCRS.GeographicCRS.Spheroid;
end

if mod(angToRotate,90)==0
    [SForward,SBackward] = internal90(angToRotate,lat,lon,Z,...
        refEllipsoid,useParallel);
else
    tmpS = internalHorizonPts(angToRotate,lat,lon,Z,...
        refEllipsoid,useParallel,method);
    [SForward,SBackward] = cleanupRotated(angToRotate,tmpS,R,method);
end
if verbose
    fprintf('rotation %f, azimuths %f %f\n',angToRotate,SForward.azm,SBackward.azm);
end
end