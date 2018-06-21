function [x,t]=serialize(s)
%SERIALIZE Serialize DBAT struct to vector.
%
%   X=SERIALIZE(S) extracts and reorders IO, EO, and OP parameters
%   from the DBAT struct S into the vector X. The parameter
%   selection and ordering is given by the S.serial field.
%
%   [X,T]=... also returns a same-sized cell array T with parameter
%   type strings.
%
%See also: DESERIALIZE, BUILDSERIALINDICES.

% Pre-allocate vector.
x=nan(s.serial.n,1);
% Insert element according to pre-calculated selection and ordering.
x(s.serial.IO.dest)=s.IO(s.serial.IO.src);
x(s.serial.EO.dest)=s.EO(s.serial.EO.src);
x(s.serial.OP.dest)=s.OP(s.serial.OP.src);

if nargout>1
    t=cell(size(x));
    t(s.serial.IO.dest)=s.paramTypes.IO(s.serial.IO.src);
    t(s.serial.EO.dest)=s.paramTypes.EO(s.serial.EO.src);
    t(s.serial.OP.dest)=s.paramTypes.OP(s.serial.OP.src);
end
