function ncwriteSpatialRef(filename,location,rasterref)
% ncwriteSpatialRef(filename,location,rasterref)
%
% filename - Net CDF file
% location - e.g. variable appropriate to this spatial reference, or '/'
% rasterref - MATLAB spatial reference object

% relevant fields of the input object
assert(contains(class(rasterref),'rasterref'),...
    'rasterref is class %s, not recognized',class(rasterref))
% world file matrix
W = worldFileMatrix(rasterref);
fn = fieldnames(rasterref);
t = contains(fn,'CRS');
wkt = wktstring(rasterref.(fn{t}));
if contains(class(rasterref),'geographic','IgnoreCase',true)
    ncwriteatt(filename,location,'Class','Geographic')
else
    ncwriteatt(filename,location,'Class','Map')
end
ncwriteatt(filename,location,'worldFileMatrix',W(:))
ncwriteatt(filename,location,'worldFileMatrixSize',size(W))
ncwriteatt(filename,location,'RasterSize',rasterref.RasterSize)
if contains(class(rasterref),'geographic','IgnoreCase',true)
    ncwriteatt(filename,location,'GeographicCRSwkt',wkt);
else
    ncwriteatt(filename,location,'ProjectedCRSwkt',wkt);
end
end