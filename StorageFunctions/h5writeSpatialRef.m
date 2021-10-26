function h5writeSpatialRef(filename,location,rasterref)
% h5writeSpatialRef(filename,location,rasterref)
%
% filename - HDF 5 file
% location - e.g. Group appropriate to this spatial reference
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
    h5writeatt(filename,location,'Class','Geographic')
else
    h5writeatt(filename,location,'Class','Map')
end
h5writeatt(filename,location,'worldFileMatrix',W)
h5writeatt(filename,location,'RasterSize',rasterref.RasterSize)
if contains(class(rasterref),'geographic','IgnoreCase',true)
    h5writeatt(filename,location,'GeographicCRSwkt',wkt{1});
else
    h5writeatt(filename,location,'ProjectedCRSwkt',wkt{1});
end
end