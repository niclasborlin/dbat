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

if nargin<2, wantedOrder={'IO','EO','OP'}; end

% Serialize each block. All x-related indices are 1-based.
[IOserial,IOdeserial,warn]=serializeblock(s.IOblock,s.estIO,s.useIOobs);
if warn
    warning(['All IO parameters in a block should be estimated or ' ...
             'fixed, not a combination. Fixed parameters will be ' ...
             'overwritten.']);
end
[EOserial,EOdeserial,warn]=serializeblock(s.EOblock,s.estEO,s.useEOobs);
if warn
    warning(['All EO parameters in a block should be estimated or ' ...
             'fixed, not a combination. Fixed parameters will be ' ...
             'overwritten.']);
end
[OPserial,OPdeserial,warn]=serializeblock(repmat(1:size(s.OP,2),3,1),s.estOP,...
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

% Store indices.
s.serial=struct('IO',IOserial,'EO',EOserial,'OP',OPserial,'n',n);
s.deserial=struct('IO',IOdeserial,'EO',EOdeserial,'OP',OPdeserial,'n',n);


% Comput serialize indices for one block
function [serial,deserial,warn]=serializeblock(block,est,useObs)

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
  [~,ia,ib]=unique([0,block(i,:)]);
  lead(i,ia(2:end)-1)=1;
end

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
