function ss=buildserialindices(ss,wantedOrder)
%BUILDSERIALINDICES Compute serialize and deserialize indices.
%
%   S0=BUILDSERIALINDICES(S0) populates the serial/deserial and
%   post.res.ix fields of the DBAT struct S0 based on the information in
%   the IO.struct.block, EO.struct.block, and bundle.est.* fields, as
%   described below. The order of the elements are IO, EO, and OP as
%   default. Use SERIALIZE(S0,ORDER), where ORDER is a cell array with
%   strings 'IO', 'EO', 'OP' to modify the order.
%
%   S0.IO.struct.block indicates how IO parameters are blocked, i.e.,
%   shared between images (columns). S0.bundle.est.IO indicates what
%   IO parameters should be estimated. The field bundle.serial.IO.src
%   will indicate what IO elements should be inserted into the vector
%   x of unknowns. The field bundle.serial.IO.dest will indicate where
%   in x the elements should be put. Only the first element from each
%   block will be copied.
%
%   Correspondingly, the bundle.deserial.IO.src indicate what
%   estimated elements in x should be copied and
%   bundle.deserial.IO.dest indicate where in IO the elements should
%   be put.
%
%   Similar reason applies to EO and OP, except the OP parameters
%   are all assumed to be distinct.
%
%   The post.res.ix struct contains subfields IP, IO, EO, OP with
%   indices into the bundle residual vector.

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
%   post.res.ix.IP - index into the residual vector for image observations.
%   post.res.ix.IO - index into the residual vector for IO observations.
%   post.res.ix.EO - index into the residual vector for EO observations.
%   post.res.ix.OP - index into the residual vector for OP observations.
%   post.res.ix.n  - total number of observations.

if nargin<2, wantedOrder={'IO','EO','OP'}; end

% Find unique IOblock and EOblock columns and indicate if columns
% are simple.
ss.IO.struct.uniq=false(1,size(ss.IO.struct.block,2));
[~,ia,ic]=unique(ss.IO.struct.block','rows');
ss.IO.struct.uniq(ia)=true;
ss.IO.struct.no=ic';
ss.IO.struct.simple=all(ss.IO.struct.block==ss.IO.struct.block(ones(end,1),:),1);

ss.EO.struct.uniq=false(1,size(ss.EO.struct.block,2));
[~,ia,ic]=unique(ss.EO.struct.block','rows');
ss.EO.struct.uniq(ia)=true;
ss.EO.struct.no=ic';
ss.EO.struct.simple=all(ss.EO.struct.block==ss.EO.struct.block(ones(end,1),:),1);

% Serialize each block. All x-related indices are 1-based.
[IOlead,IOserial,IOdeserial,warn,blockIx]=serializeblock(ss.IO.struct.block,...
                                                  ss.bundle.est.IO, ...
                                                  ss.IO.prior.use);
if warn
    warning(['All IO parameters in a block should be estimated or ' ...
             'fixed, not a combination. Fixed parameters will be ' ...
             'overwritten.']);
end
switch length(blockIx)
  case 0
    % No IO blocks. Legacy models will work with one IO column per
    % image. FIXME
    ss.EO.cam=1:size(ss.EO.val,2);
  case 1
    % One IO block. Legacy models require column indices to the
    % block. WARNING: Untested for blockIx!=1.
    ss.EO.cam=repmat(blockIx,1,size(ss.EO.val,2));
  otherwise
    % More than one code block => signal incompatibility with NaN's.
    ss.EO.cam=nan(1,size(ss.EO.val,2));
end    
[EOlead,EOserial,EOdeserial,warn]=serializeblock(ss.EO.struct.block,...
                                                 ss.bundle.est.EO,...
                                                 ss.EO.prior.use);
if warn
    warning(['All EO parameters in a block should be estimated or ' ...
             'fixed, not a combination. Fixed parameters will be ' ...
             'overwritten.']);
end
[~,OPserial,OPdeserial,warn]=serializeblock(repmat(1:size(ss.OP.val,2),3,1),...
                                            ss.bundle.est.OP,...
                                            ss.OP.prior.use);
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

% Store arrays indicating leading parameters, i.e., the first
% parameter within each block to be estimated.
ss.IO.struct.lead=IOlead;
ss.EO.struct.lead=EOlead;

% Mask out any useObs that correspond to a repeated parameter.
ss.IO.prior.use=ss.IO.prior.use & ss.IO.struct.lead;
ss.EO.prior.use=ss.EO.prior.use & ss.EO.struct.lead;

% Store indices.
ss.bundle.serial=struct('IO',IOserial,...
                        'EO',EOserial,...
                        'OP',OPserial,...
                        'n',n);
ss.bundle.deserial=struct('IO',IOdeserial,...
                          'EO',EOdeserial,...
                          'OP',OPdeserial,...
                          'n',n);

% Indices for observations.
% Create indices into the residual vector. nObs is the total number
% of observations.
numObs=[nnz(ss.IP.vis)*2,length(IOserial.obs),length(EOserial.obs),...
        length(OPserial.obs)];
[obsIPix,obsIOix,obsEOix,obsOPix,nObs]=indvec(numObs);

ss.post.res.ix=struct('IP',obsIPix,...
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
