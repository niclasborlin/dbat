function s=deserialize(s,x,i,what)
%DESERIALIZE Update DBAT struct from vector.
%
%   S=DESERIALIZE(S,X) updates IO, EO, and OP parameters in the DBAT
%   struct S from the vector X. The parameter selection and ordering
%   is given by the S.deserial field.
%
%   S=DESERIALIZE(S,E,I) takes the updates from iteration I in the
%   bundle iteration struct E instead. I=0 corresponds to the initial
%   values. The Inf value can be used to get the final iteration.
%
%   IO=DESERIALIZE(S,E,V,'IO') returns the IO values only for the
%   (zero-based) iteration numbers stored in V. Use V='all' to get all
%   iterations. The returned IO has the same size as s.IO but extended
%   in the third dimension. The returned IO(:,:,i) values will be the
%   IO values at iteration V(i).
%
%   EO=DESERIALIZE(S,E,V,'EO') returns the EO values instead.
%   OP=DESERIALIZE(S,E,V,'OP') returns the OP values instead.
%
%See also: SERIALIZE, BUILDSERIALINDICES.

switch nargin
  case 2
    % DESERIALIZE(S,X).
    
    % Copy elements to pre-calculated selection and ordering.
    s.IO.val(s.bundle.deserial.IO.dest)=x(s.bundle.deserial.IO.src);
    s.EO.val(s.bundle.deserial.EO.dest)=x(s.bundle.deserial.EO.src);
    s.OP.val(s.bundle.deserial.OP.dest)=x(s.bundle.deserial.OP.src);
  case 3
    % DESERIALIZE(S,E,I).
    E=x;

    if isempty(E)
        error('Empty bundle iteration struct');
    end
    
    if i==inf
        % Keep I zero-based.
        i=size(E.trace,2)-1;
    end
    % Copy elements to pre-calculated selection and ordering.
    s.IO.val(s.bundle.deserial.IO.dest)=E.trace(s.bundle.deserial.IO.src,i+1);
    s.EO.val(s.bundle.deserial.EO.dest)=E.trace(s.bundle.deserial.EO.src,i+1);
    s.OP.val(s.bundle.deserial.OP.dest)=E.trace(s.bundle.deserial.OP.src,i+1);
  case 4
    % DESERIALIZE(S,E,V,'IO'/'EO'/'OP')
    E=x;
    v=i;
    
    if ischar(v) && strcmp(v,'all')
        if isempty(E)
            v=zeros(1,0);
        else
            v=0:size(E.trace,2)-1;
        end
    end
    
    if isempty(E) && ~isempty(v)
        error('Empty bundle iteration struct');
    end
    
    switch what
      case 'IO'
        % Replicate IO values.
        IO=repmat(s.IO.val,[1,1,length(v)]);
        if ~isempty(v)
            destIx=s.bundle.deserial.IO.dest;
            srcIx=s.bundle.deserial.IO.src;
            % Convert destination indices to row, column.
            [ii,jj]=ind2sub(size(s.IO.val),destIx);
            for i=1:length(v)
                % Update layer i.
                IO(sub2ind(size(IO),ii,jj,repmat(i,size(ii))))=...
                    E.trace(srcIx,v(i)+1);
            end
        end
        s=IO;
      case 'EO'
        % Replicate EO values.
        EO=repmat(s.EO.val,[1,1,length(v)]);
        if ~isempty(v)
            destIx=s.bundle.deserial.EO.dest;
            srcIx=s.bundle.deserial.EO.src;
            % Convert destination indices to row, column.
            [ii,jj]=ind2sub(size(s.EO.val),destIx);
            for i=1:length(v)
                % Update layer i.
                EO(sub2ind(size(EO),ii,jj,repmat(i,size(ii))))=...
                    E.trace(srcIx,v(i)+1);
            end
        end
        s=EO;
      case 'OP'
        % Replicate OP values.
        OP=repmat(s.OP.val,[1,1,length(v)]);
        if ~isempty(v)
            destIx=s.bundle.deserial.OP.dest;
            srcIx=s.bundle.deserial.OP.src;
            % Convert destination indices to row, column.
            [ii,jj]=ind2sub(size(s.OP.val),destIx);
            for i=1:length(v)
                % Update layer i.
                OP(sub2ind(size(OP),ii,jj,repmat(i,size(ii))))=...
                    E.trace(srcIx,v(i)+1);
            end
        end
        s=OP;
      otherwise
        error('Bad selector');
    end
end
