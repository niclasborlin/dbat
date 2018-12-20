function s=parseblockvariant(s)
%PARSEBLOCKVARIANT Parse IO/EO block-variant info of a DBAT structure.
%
%   S=PARSEBLOCKVARIANT(S) parses the block-variant IO and EO
%   information in S.

%   Parse the S.IO.struct.block and S.EO.struct.block fields are
%   parsed and populate the uniq, no, and isSimple fields. The
%   leading field is populated by buildserialindices.

% Find unique IOblock and EOblock columns and indicate if columns
% are simple.
[s.IO.struct.uniq,s.IO.struct.no,s.IO.struct.isSimple]=...
    parseblock(s.IO.struct.block);

[s.EO.struct.uniq,s.EO.struct.no,s.EO.struct.isSimple]=...
    parseblock(s.EO.struct.block);


function [uniq,no,isSimple]=parseblock(block)

% Find out which columns are unique.
uniq=false(1,size(block,2));
[~,ia,ic]=unique(block','rows');
uniq(ia)=true;
% Give each column a number among the unique
no=ic';
isSimple=all(block==block(ones(end,1),:),1);
