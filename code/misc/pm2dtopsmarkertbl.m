function psTbl=pm2dtopsmarkertbl(pmTbl,varargin)
%PM2DTOPSMARKERTBL Convert Photomodeler 2D table to Photoscan marker table.
%
%   PSPTS=PM2DTOPSMARKERTBL(PMPTS) converts the Photomodeler 2D
%   table PMPTS to a Photoscan marker table with fields
%       markerId - N-array with marker ids copied from PMPTS.id.
%       cameraId - N-array with camera ids copied as PMPTS.imNo-1.
%       pos      - 2-by-N array with image coordinates copied from PMPTS.pos.
%       pinned   - logical N-array with pinned status. Defaults to true.
%
%   PSPTS=PM2DTOPSMARKERTBL(PMPTS,MAP) uses the sparse array MAP to
%   map Photomodeler ids to Photoscan.
%
%   ...=PM2DTOPSMARKERTBL(...,B) uses the bool scalar or vector B
%   to set the pinned status.


map=[]
b=true;

for i=1:length(varargin)
    if islogical(varargin{i})
        b=varargin{i};
    else
        map=varargin{i};
    end
end

% Scalar expansion of b.
if isscalar(b), b=repmat(b,1,length(pmTbl.id)); end

if isempty(map)
    markerId=pmTbl.id;
else
    markerId=full(map(pmTbl.id));
end

cameraId=pmTbl.imNo-1;

pos=pmTbl.pos;

pinned=b;

psTbl=struct('markerId',markerId,'cameraId',cameraId,'pos',pos,'pinned',pinned);