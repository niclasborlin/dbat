function s=buildserialindices(s,wantedOrder)
%BUILDSERIALINDICES Compute serialize and deserialize indices.
%
%   S0=BUILDSERIALINDICES(S0) populates the serial/deserial fields of
%   the DBAT struct S0 based on the information in the IOblock,
%   EOblock, and estIO, estEO, estOP fields, as described below. The
%   order of the elements are IO, EO, and OP as default. Use
%   SERIALIZE(S0,ORDER), where ORDER is a cell array with strings
%   'IO', 'EO', 'OP' to modify the order.
%
%   S0.IOblock indicates how IO parameters are blocked, i.e., shared
%   between images (columns). S0.estIO indicates what IO parameters
%   should be estimated. The field serial.IO.src will indicate what IO
%   elements should be inserted into the vector x of unknowns. The
%   field serial.IO.dest will indicate where in x the elements should
%   be put. Only the first element from each block will be copied.
%
%   Correspondingly, the deserial.IO.src indicate what estimated
%   elements in x should be copied and deserial.IO.dest indicate where
%   in IO the elements should be put.
%
%   Similar reason applies to EO and OP, except the OP parameters
%   are all assumed to be distinct.

% Serialization:
%
%   serial.IO.src is index into IO.
%   serial.IO.dest is index into x.
%   serial.IO.obs is index into serial.IO.src and serial.IO.dest.
%   serial.EO.src is index into EO.
%   serial.EO.dest is index into x.
%   serial.EO.obs is index into serial.EO.src and serial.EO.dest.
%   serial.OP.src is index into OP.
%   serial.OP.dest is index into x.
%   serial.OP.obs is index into serial.OP.src and serial.OP.dest.
%
% Deserialization:
%   deserial.IO.src is index into x.
%   deserial.IO.dest is index into IO.
%   deserial.EO.src is index into x.
%   deserial.EO.dest is index into EO.
%   deserial.OP.src is index into x.
%   deserial.OP.dest is index into OP.
%
% Observation ordering
%   residuals.ix.IP - index into the residual vector for image observations.
%   residuals.ix.IO - index into the residual vector for IO observations.
%   residuals.ix.EO - index into the residual vector for EO observations.
%   residuals.ix.OP - index into the residual vector for OP observations.
%   residuals.ix.n  - total number of observations.

if nargin<2, wantedOrder={'IO','EO','OP'}; end

% Find unique IOblock and EOblock columns and indicate if columns
% are simple.
s.IOunique=false(1,size(s.IOblock,2));
[~,ia]=unique(s.IOblock','rows');
s.IOunique(ia)=true;
s.IOsimple=all(s.IOblock==s.IOblock(ones(end,1),:),1);

s.EOunique=false(1,size(s.EOblock,2));
[~,ia]=unique(s.EOblock','rows');
s.EOunique(ia)=true;
s.EOsimple=all(s.EOblock==s.EOblock(ones(end,1),:),1);

% Serialize each block. All x-related indices are 1-based.
[IOlead,IOserial,IOdeserial,warn,blockIx]=serializeblock(s.IOblock,s.estIO,...
                                                  s.useIOobs);
if warn
    warning(['All IO parameters in a block should be estimated or ' ...
             'fixed, not a combination. Fixed parameters will be ' ...
             'overwritten.']);
end
switch length(blockIx)
  case 0
    % No IO blocks. Legacy models will work with one IO column per
    % image.
    s.imCams=1:size(s.EO,2);
  case 1
    % One IO block. Legacy models require column indices to the
    % block. WARNING: Untested for blockIx!=1.
    s.imCams=repmat(blockIx,1,size(s.EO,2));
  otherwise
    % More than one code block => signal incompatibility with NaN's.
    s.imCams=nan(1,size(s.EO,2));
end    
[EOlead,EOserial,EOdeserial,warn]=serializeblock(s.EOblock,s.estEO,s.useEOobs);
if warn
    warning(['All EO parameters in a block should be estimated or ' ...
             'fixed, not a combination. Fixed parameters will be ' ...
             'overwritten.']);
end
[~,OPserial,OPdeserial,warn]=serializeblock(repmat(1:size(s.OP,2),3,1),s.estOP,...
                                          s.useOPobs);
if warn
    warning(['All OP parameters in a block should be estimated or ' ...
             'fixed, not a combination. Fixed parameters will be ' ...
             'overwritten.']);
end

% Adjust indices to correspond to different parts of the x vector.
n=0;
for i=1:length(wantedOrder)
    switch wantedOrder{i}
      case 'IO'
        % Adjust references into x.
        IOserial.dest=IOserial.dest+n;
        IOdeserial.src=IOdeserial.src+n;
        n=n+length(IOserial.dest);
      case 'EO'
        % Adjust references into x.
        EOserial.dest=EOserial.dest+n;
        EOdeserial.src=EOdeserial.src+n;
        n=n+length(EOserial.dest);
      case 'OP'
        % Adjust references into x.
        OPserial.dest=OPserial.dest+n;
        OPdeserial.src=OPdeserial.src+n;
        n=n+length(OPserial.dest);
    end
end

% Storead leading arrays.
s.IOlead=IOlead;
s.EOlead=EOlead;

% Store indices.
s.serial=struct('IO',IOserial,'EO',EOserial,'OP',OPserial,'n',n);
s.deserial=struct('IO',IOdeserial,'EO',EOdeserial,'OP',OPdeserial,'n',n);

% Indices for observations.
% Create indices into the residual vector. nObs is the total number
% of observations.
numObs=[nnz(s.vis)*2,length(IOserial.obs),length(EOserial.obs),...
        length(OPserial.obs)];
[obsIPix,obsIOix,obsEOix,obsOPix,nObs]=indvec(numObs);

s.residuals.ix=struct('IP',obsIPix,...
                      'IO',obsIOix,...
                      'EO',obsEOix,...
                      'OP',obsOPix,...
                      'n',nObs);

% Compute serialize indices for one block
function [lead,serial,deserial,warn,blockIx]=serializeblock(block,est,useObs)

% Elements within one parameter block should all be marked as
% 'estimate' or 'fixed', not a combination.
warn=false;
for i=1:size(block,1)
    estBlock=unique(block(i,est(i,:)));
    fixedBlock=unique(block(i,~est(i,:)));
    if ~isempty(intersect(estBlock,fixedBlock))
        warn=true;
    end
end

% Do not touch fixed elements.
block(~est)=0;

% Find the leading elements of each row.
lead=zeros(size(block));
for i=1:size(block,1)
  [~,ia,~]=unique([0,block(i,:)]);
  lead(i,ia(2:end)-1)=1;
end

% Return columns corresponding to parameter blocks.
blockIx=find(any(lead,1));

% Indices to serialize v (matrix -> vector)
serial.src=find(lead);
serial.dest=(1:nnz(lead))';

% Find subset of serial indices that correspond to parameters with
% prior observations.
serial.obs=find(useObs(lead>0));

% Now create the inverse mapping to distribute x elements to IO.
% (This can be done quicker.)
dist=lead;
dist(lead~=0)=serial.dest;

% Expand to all elements that should be equal.
for k=1:length(serial.dest)
    [i,j]=find(dist==k);
    % Distribute over the block.
    inBlock=block(i,:)==block(i,j);
    dist(i,inBlock)=k;
end

% Extract indices.
[i,~,v]=find(dist(:));
deserial.dest=i;
deserial.src=v;
