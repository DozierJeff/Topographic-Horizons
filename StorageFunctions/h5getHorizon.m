function [ H, varargout ] = h5getHorizon( h5file, varargin )
% angles = h5getHorizon(h5file)
% angles = h5getHorizon(h5file,azimuth)
%
% retrieve horizon angles or horizon mask
%
% Input
%   h5file - HDF 5 filename
% Optional input, in order
%   azimuth - in degrees, either scalar or matrix of size of topographic grid,
%       for the direction(s) in which the horizon angles are desired
%       (orientation set by azimuthPreference, either azimuth to south is 0
%       degrees, + east of south, - west of south, or azimuth clockwise
%       from north
%
% Output
%   if just azimuth input argument, angles to the horizon in that direction
%   otherwise, horizon angles in all directions
% Optional output
%   azimuth(s) corresponding to horizons

p = inputParser;
addRequired(p,'h5file',@ischar)
addOptional(p,'azimuth',[],@isnumeric)
parse(p,h5file,varargin{:});
azm = p.Results.azimuth;

% check to make sure input is a .h5 file
[~,~,ext] = fileparts(h5file);
assert(strcmpi(ext,'.h5'),'input must be HDF 5 file')
info = h5info(h5file);
siz = info.Groups(1).Datasets(3).ChunkSize;
FillValue = info.Groups(1).Datasets(3).FillValue;
if ismatrix(siz)
    siz = [siz(1) siz(2) 1];
end

% heading information
angles = h5readatt(h5file,'/Grid','azimuths');

% horizon angles, get all if input azimuth empty
if isempty(azm)
    H = h5read(h5file,'/Grid/horizons');
elseif isscalar(azm)
    adiff = abs(angles-azm);
    idx = find(adiff==min(adiff));
    H = h5read(h5file,'/Grid/horizons',[1 1 idx],siz);
else
    assert(isequal(size(azm),siz(1:2)),...
        'if azimuth is not a scalar, it must be size of horizon grid')
    u = unique(azm(:));
    for k=1:length(u)
        adiff = abs(angles-u(k));
        idx = find(adiff==min(adiff));
        H = h5read(h5file,'/Grid/horizons', [1 1 idx], siz);
        if k==1
            X = zeros(size(H),'like',H);
        end
        t = azm==u(k);
        X(t) = H(t);
    end
    % all of output should be filled, except for cells where sun has not
    % yet risen
    if ~isnan(FillValue)
    X(output==FillValue) = NaN;
    end
    H = X;
end
if nargout>1
    if isscalar(azm)
        varargout{1} = angles(idx);
    else
        varargout{1} = angles;
    end
end