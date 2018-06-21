function s=deserialize(s,x)
%DESERIALIZE Update DBAT struct from vector.
%
%   S=DESERIALIZE(S,X) updates IO, EO, and OP parameters in the DBAT
%   struct S from the vector X. The parameter selection and ordering
%   is given by the S.deserial field.
%
%See also: SERIALIZE, BUILDSERIALINDICES.

% Copy elements to pre-calculated selection and ordering.
s.IO(s.deserial.IO.dest)=x(s.deserial.IO.src);
s.EO(s.deserial.EO.dest)=x(s.deserial.EO.src);
s.OP(s.deserial.OP.dest)=x(s.deserial.OP.src);
