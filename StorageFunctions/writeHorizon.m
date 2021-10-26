function [outputFile] = writeHorizon(filename,R,azm,horizons,varargin)
% outputFile(s) = writeHorizon(filename,R,azm,horizons [,distances])
%save horizons in format specified by filename:
%options are geotiff, HDF 5, NetCDF 4, or MATLAB
%if geotiff, only option is to read the whole file but then can be
%   processed by other software
%if HDF 5 (.h5) or NetCDF 4, output is stored in block compressed form so
%   horizons for a specific azimuth can be retrieved without reading the whole file
%if MATLAB (.mat), output is stored as an interpolating function so that
%   horizons can be retrieved for specific azimuths and locations
%   (to wrute just the azimuths and horizons to a MATLAB file, use the
%   MATLAB save function)
%Note: this routine works only for R2020b or later, because it depends on
%the ProjectedCRS or GeographicCRS fields in R to be populated. Sorry.
%Please contact me (dozier@ucsb.edu) if this limitation affects you,
%because I can code an alternative.
%
%input
% filename, either .tif, .h5, .nc, or .mat, include full path if needed (if
%   output file is a .tif and distance output is specified, two output files
%   will be created, one filename containing 'horizon' and the other 'distance')
% R - raster reference for the horizon grid, can be a MapCellsReference
%   or MapPostingsReference object, or a GeographicCellsReference or
%   GeographicPostingsReference object
%   (If it's a 'postings' object, will be converted to 'cells')
% azm - vector of azimuths of horizons in degrees either from -180 to +180
%   ccw from north, with 0 south, or cw from 0 to 360 with 0 north, depending
%   on how the azimuthPreference function is set (values from horizonAllDirections)
% horizons - horizons in directions azm as BSQ array (output from
%   horizonAllDirections), one x-y or lat-lon plane for each azimuth, as
%   floating point numbers
%optional input, name-value pairs
% 'dis' - distances in directions azm in BSQ format (output from
%   horizonAllDirections), one x-y or lat-lon plane for each azimuth, as
%   floating point numbers, must be same size as horizons array
% 'GeoKeyDirectoryTag' - only if output file(s) in .tif format, and only if
%   the program cannot figure out the projection or geographic parameters
%   from the input raster reference object

assert(~verLessThan('matlab','9.9'),...
    'this function works only on MATLAB R2020b or later')

p = inputParser;
addRequired(p,'filename',@ischar)
addRequired(p,'R',@(x) contains(class(x),'rasterref'))
addRequired(p,'azm',@(x) isnumeric(x) && isvector(x))
addRequired(p,'horizons',@(x) isnumeric(x) && ndims(x)==3)
addParameter(p,'dis',[],@(x) isnumeric(x) && ndims(x)==3)
addParameter(p,'GeoKeyDirectoryTag',struct([]),@isstruct)
parse(p,filename,R,azm,horizons,varargin{:})

% file must be .tif, .h5, .nc, or .mat
[~,~,ext] = fileparts(filename);
geotiff = false;
hdf5 = false;
mat = false;
netCDF = false;
switch ext
    case '.tif'
        geotiff = true;
    case '.h5'
        hdf5 = true;
    case '.mat'
        mat = true;
    case '.nc'
        netCDF = true;
    otherwise
        error('filename extension must be .tif, .h5, .nc, or .mat')
end
% distances must be same size if entered
if ~isempty(p.Results.dis)
    distances = p.Results.dis;
    assert(isequal(size(distances),size(horizons)),...
        'if specified, ''dis'' must be same size as horizons')
end
% change R to 'cells' if currently 'postings
R = postings2cells(R);

assert(isfloat(horizons),...
    'input horizon array must be floating point (I could revise the code)')
horizons = single(horizons);
if exist('distances','var')
    distances = single(distances);
end

if hdf5
    outputFile = writeHDF5;
elseif geotiff
    outputFile = writegeotiff;
elseif mat
    outputFile = writeMAT;
elseif netCDF
    outputFile = writeNetCDF;
end

    function matname=writeMAT()
        % file name should include 'horizon'
        [fpath,fname,xt] = fileparts(filename);
        fullname = fullfile(fpath,[fname xt]);
        %create lookup function
        [r,c,a] = ndgrid(1:size(horizons,1),1:size(horizons,2),azm);
        Fhorizon = griddedInterpolant(r,c,a,horizons,'linear','nearest');
        % same for distances
        if exist('distances','var')
            % interpolating function
            Fdistance = griddedInterpolant(r,c,a,distances,'linear','nearest');
        end
        % write file
        if exist('Fdistance','var')
            save(fullname,'-v7.3','R','Fhorizon','Fdistance')
        else
            save(fullname,'-v7.3','R','Fhorizon')
        end
        matname = fullname;
    end

    function gname = writegeotiff(varargin)
        pg = inputParser;
        addOptional(pg,'GeoKeyDirectoryTag',struct([]),@isstruct)
        parse(pg,varargin{:})
        [fpath,fname,xt] = fileparts(filename);
        fullname = fullfile(fpath,[fname xt]);
        if isempty(pg.GeoKeyDirectoryTag)
            if contains(class(R),'geographic','IgnoreCase',true)
                geotiffwrite(fullname,horizons,R)
            else
                try
                    wkt = char(wktstring(R.ProjectedCRS));
                    e = strfind(wkt,'EPSG');
                    lastpart = wkt(e(end):end);
                    kg = strfind(lastpart,',');
                    lastpart = lastpart(kg+1:end);
                    kg = strfind(lastpart,']');
                    code = str2double(lastpart(1:kg-1));
                    geotiffwrite(fullname,horizons,R,'CoordRefSysCode',code)
                catch
                    error('check input, you may need to include ''GeoKeyDirectoryTag''')
                end
            end
        else
            geotiffwrite(fullname,horizons,R,...
                'GeoKeyDirectoryTag',pg.GeoKeyDirectoryTag)
        end
        gname = fullname;

        % same for distances, add 'distance' to filename
        if exist('distances','var')
            fname = [fname '_distance'];
            disname = fullfile(fpath,[fname xt]);
            if isempty(pg.GeoKeyDirectoryTag)
                if contains(class(R),'geographic','IgnoreCase',true)
                    geotiffwrite(disname,distances,R)
                else
                    geotiffwrite(disname,distances,R,'CoordRefSysCode',code)
                end
            else
                geotiffwrite(disname,distances,R,...
                    'GeoKeyDirectoryTag',pg.GeoKeyDirectoryTag)
            end
            gname = {gname,disname}.';
        end
    end

    function h5file = writeHDF5()
        % horizons HDF specific stuff
        [fpath,fname,xt] = fileparts(filename);
        h5file = fullfile(fpath,[fname xt]);
        deflateLevel = 9; % compression factor
        group = '/Grid';
        fillValue = nan(1,1,'like',horizons);

        % create file and write data
        dname = 'horizons';
        h5create(h5file,[group '/' dname],size(horizons),...
            'Deflate',deflateLevel,...
            'ChunkSize',[size(horizons,1) size(horizons,2) 1],...
            'FillValue',fillValue,...
            'DataType',class(horizons));
        h5write(h5file,[ group '/' dname],horizons)
        h5writeatt(h5file,[ group '/' dname],'Units','degree');
        if exist('distances','var')
            dname = 'distances';
            fillValue = nan(1,1,'like',distances);
            h5create(h5file,[group '/' dname],size(distances),...
                'Deflate',deflateLevel,...
                'ChunkSize',[size(distances,1) size(distances,2) 1],...
                'FillValue',fillValue,...
                'DataType',class(distances));
            h5write(h5file,[ group '/' dname],distances)
            h5writeatt(h5file,[ group '/' dname],'Units','meter');
        end

        % write Global Attributes
        h5writeatt(h5file,group,'azimuths',azm(:).');
        h5create(h5file,[group '/' 'SpatialRef'],GetSize(R))
        h5writeSpatialRef(h5file,[group '/' 'SpatialRef'],R)
    end

    function ncfile=writeNetCDF()
        % horizons NetCDF specific stuff
        [fpath,fname,xt] = fileparts(filename);
        ncfile = fullfile(fpath,[fname xt]);
        deflateLevel = 9; % compression factor
        fillValue = nan(1,1,'like',horizons);

        % create file and write data
        group = '/Grid';
        dname = 'horizons';
        nccreate(ncfile,[group '/' dname],...
            'Deflate',deflateLevel,...
            'ChunkSize',[size(horizons,1) size(horizons,2) 1],...
            'FillValue',fillValue,...
            'DataType',class(horizons),...
            'Dimensions',{'row',size(horizons,1),'col',size(horizons,2),...
            'azm',length(azm)});
        ncwrite(ncfile,[group '/' dname],horizons)
        ncwriteatt(ncfile,[group '/' dname],'Units','degree')
        if exist('distances','var')
            dname = 'distances';
            nccreate(ncfile,[group '/' dname],...
                'Deflate',deflateLevel,...
                'ChunkSize',[size(distances,1) size(distances,2) 1],...
                'FillValue',fillValue,...
                'DataType',class(distances),...
                'Dimensions',{'row',size(distances,1),'col',size(distances,2),...
                'azm',length(azm)});
            ncwrite(ncfile,[group '/' dname],distances)
            ncwriteatt(ncfile,[group '/' dname],'Units','meter')
        end

        % write Global Attributes
        ncwriteatt(ncfile,'/','azimuths',azm(:).')
        ncwriteSpatialRef(ncfile,'/',R)
    end
end