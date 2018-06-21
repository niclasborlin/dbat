function s=serialize(s,order)
%SERIALIZE Compute serialize and deserialize indices.
%
%   S0=SERIALIZE(S0) populates the serial/deserial fields of the DBAT
%   struct S0 based on the information in the IOblock, EOblock, and
%   estIO, estEO, estOP fields, as described below. The order of the
%   elements are IO, EO, and OP as default. Use SERIALIZE(S0,ORDER),
%   where ORDER is a cell array with strings 'IO', 'EO', 'OP' to
%   modify the order.
%
%   S0.IOblock indicates how IO parameters are blocked, i.e., shared
%   between images (columns). S0.estIO indicates what IO parameters
%   should be estimated. The field serial.IOIO will indicate what IO
%   elements should be inserted into the vector x of unknowns. The
%   field serial.IOx will indicate where in x the elements should be
%   put. Only the first element from each block will be copied.
%
%   Correspondingly, the deserial.IOx indicate what estimated elements
%   in x should be copied and deserial.IOIO indicate where in IO
%   the elements should be put.
%
%   Similar reason applies to EO and OP, except the OP parameters
%   are all assumed to be distinct.

if nargin<2, order={'IO','EO','OP'}; end

% Extract what parameters of IOblock to estimate.
block=s.IOblock;
% Do not touch fixed elements.
block(~s.estIO)=0;

% Find the leading elements of each row.
lead=zeros(size(block));
for i=1:size(block,1)
  [~,ia,ib]=unique([0,block(i,:)]);
  lead(i,ia(2:end)-1)=true;
end

% Indeces to serialize v (matrix -> vector)
serial.IOIO=find(lead);
serial.IOx=1:nnz(lead);

% Now create the inverse mapping to distribute x elements to IO.
% (This can be done quicker.)
dist=lead;
dist(lead~=0)=serial.IOx;

% Expand to all elements that should be equal.
for k=1:length(serial.IOx)
    [i,j]=find(lead==k);
    % Distribute over the block.
    inBlock=block(i,:)==block(i,j);
    dist(i,inBlock)=k;
end

[i,j,v]=find(dist(:));
deserial.IOIO=i;
deserial.IOx=v;
