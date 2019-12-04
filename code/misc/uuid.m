function id=uuid
%UUID Return a Universal Unique IDentifier.
%
%   ID=UUID calls the randomUUID() java function to compute a
%   random UUID.

id=char(java.util.UUID.randomUUID());
