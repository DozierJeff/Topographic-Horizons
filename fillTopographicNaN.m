function newZ = fillTopographicNaN(Z,varargin)
% newZ = fillTopographicNaN(Z [,method])
%
%fill topographic NaNs in DEM that has legitimate NaNs around the edges of
%the area of interest but has NaNs in the interior that need to be filled
%(might be applicable to other imagery with NaNs)
%
%The optional argument specifies the method to fill the NaNs, choices (can
%truncate to first 3 letters) are:
%   'regionfill' - default, uses regionfill
%   'coherent' - uses inpaintCoherent
%   'exemplar' - uses inpaintExemplar

p = inputParser;
addRequired(p,'Z',@(x) isnumeric(x) && (ismatrix(x) || ndims(x)==3))
addOptional(p,'method','regionfill',@(x) ischar(x) || isstring(x))
parse(p,Z,varargin{:})

Z = p.Results.Z;
method = lower(p.Results.method);

% identify the NaNs to fill (i.e., only some of the NaNs)
allNaN = isnan(Z);
boundaryNaN = ~imfill(~allNaN,8,'holes');
NaNtoFill = allNaN & ~boundaryNaN;

% fill the NaNs
Z(allNaN) = 0;
if strncmp(method,'reg',3)
    newZ = regionfill(Z,NaNtoFill);
elseif strncmp(method,'coh',3)
    newZ = inpaintCoherent(Z,NaNtoFill);
elseif strncmp(method,'exe',3)
    newZ = inpaintExemplar(Z,NaNtoFill);
else
    error('optional argument ''%s'' not recognized',p.Results.method)
end
newZ(boundaryNaN) = NaN;

end