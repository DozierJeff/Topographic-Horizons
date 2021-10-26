function [SForward,SBackward] = cleanupRotated(angleRotated,tmpS,R,method)
%cleanup results from internalHorizonPts, putting values back in original grid
%size and orientation
%
%Input
%   tmpS - output from internalHorizonPts
%   R - raster reference object of original elevation grid
%   method - of the original rotation, 'nearest' is the mechanism currently
%       supported, 'bilinear' is not but might be in the future
%
%Output
%   SForward and SBackward contain the mean azimuth and grids of horizon
%       angles and distances to the horizon in the forward and backward
%       directions
%
%Process
%   The original rotation leaves out some pixels, duplicates some, and can
%   result in a ragged edge. The routine rotates the original row and
%   column coordinates and uses their position in the rotated grid to fill
%   values in the new grid of the same size as the original un-rotated.
%   Missing values, which will consitute up to 18% for rotation of 45 deg,
%   are filled in by coherent inpainting.

% rows and columns of the orginal grid ne
assert(strcmpi(method,'nearest'),...
    'this routine works only on rotations via nearest-neighbor')
[row,col] = ndgrid(1:R.RasterSize(1),1:R.RasterSize(2));
% rotated, and made column vectors
rr = imrotate(row,angleRotated);
rc = imrotate(col,angleRotated);
rr = rr(:);
rc = rc(:);
% indexes to only values in original grid
tj = rr~=0;
idxJ = find(tj)';

SForward.azm = tmpS.Forward.azm;
SBackward.azm = tmpS.Backward.azm;

interpV = [];
substruct = {'Forward','Backward'};
vName = {'horzAng','horzDis'};
for k=1:length(substruct)
    for n=1:length(vName)
        tmpV = tmpS.(substruct{k}).(vName{n});
        tmpV = tmpV(:);
        % tmpOut is the output variable, but with some missing values
        tmpOut = nan(size(row,1),size(row,2));
        for j=idxJ
            tmpOut(rr(j),rc(j)) = tmpV(j);
        end
        % inpainting of missing values
        tfill = isnan(tmpOut);
        if any(tfill,'all')
            % fill-in mask can't be NaN, but none of the variables are negative
            % so can use a negative number to flag
            tmpOut(tfill) = -100;
            tmpOut = inpaintCoherent(tmpOut,tfill);
            % use different method if we haven't got them all (usually we have)
            tfill = tmpOut<0;
            if any(tfill,'all')
                tmpOut = regionfill(tmpOut,tfill);
            end
        end
        interpV = cat(3,interpV,tmpOut);
    end
end

% output structures
SForward.horzAng = single(interpV(:,:,1));
SForward.horzDis = single(interpV(:,:,2));
SBackward.horzAng = single(interpV(:,:,3));
SBackward.horzDis = single(interpV(:,:,4));
end