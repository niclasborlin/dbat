function ret=buildparamtypes(s,sel)
%BUILDPARAMTYPES Build parameter type strings.
%
%   S1=BUILDPARAMTYPES(S0) populates the IO.type, EO.type, and OP.type fields
%   of the DBAT struct S1 based on the information in the DBAT struct S0, as
%   described below.
%
%   T=BUILDPARAMTYPES(S,SEL), where SEL is 'IO', 'EO', or 'OP', returns the
%   cell array of selected parameter types instead of the full struct.
%
%   Each S.IO column stores the parameters below.
%   - cc      - camera constant in camera units.
%   - px, py  - principal point in camera units.
%   - as      - affine parameter. Aspect will be (1+af):1.
%   - sk      - skew parameters.
%   - K1,...  - radial distortion parameters of Brown (1971).
%   - P1,...  - tangential distortion parameters of Brown (1971).
%   The number of K and P parameter are given by s.IO.model.nK and
%   s.IO.model.nP.
%
%   The names above are returned in the NC-by-nCams IO cell array, where
%   NC=5+s.IO.model.nK+s.IO.model.nP. If multiple S.IO columns are present,
%   the column number is appended. If the block number is different from the
%   image number, the block number is appended in parantheses.
%
%   Each S.EO column stores the parameters below. The first 6 parameters
%   may be estimated by the bundle.
%   - EX, EY, EZ        - external coordinates of the camera center in
%                         project units. 
%   - omega, phi, kappa - Euler angles for the camera orientation in radians.
%
%   The first two letters of the names above are returned in the 6-by-nImages
%   EO cell array. If multiple EO columns are present, a camera identifier is
%   appended to each parameter. The camera identifier consists of the camera
%   sequence number and camera id.
%
%   Each S.OP column stores the X, Y, Z coordinates. The OP names are returned
%   in the 3-by-nOP cell array OP. Object points are prefixed with 'O', i.e.
%   'OX', 'OY', 'OZ'. Control points are prefixed with 'C'. Check points are
%   prefixed with 'H'. Furthermore, a point identifier is appended. The point
%   identifier consists of the sequence number and if necessary, the OP id,
%   the OP raw id, and the OP label.
%
%See also: PROB2DBATSTRUCT.

if nargin<2, sel='struct'; end

IOtypes={};
EOtypes={};
OPtypes={};

if strcmp(sel,'IO') || strcmp(sel,'struct')
    % IO parameter types.
    Knames=arrayfun(@(x)sprintf('K%d',x),1:s.IO.model.nK,'uniformoutput',false);
    Pnames=arrayfun(@(x)sprintf('P%d',x),1:s.IO.model.nP,'uniformoutput',false);
    IOtypes={'cc','px','py','as','sk',Knames{:},Pnames{:}}';
    if size(s.IO.val,2)>1
        IOtypes=repmat(IOtypes,1,size(s.IO.val,2));
        if isscalar(unique(s.IO.struct.block(:)))
            % All in one block, do nothing.
        elseif all(s.IO.struct.isSimple)
            % Purely image-invariant.
            for i=1:size(s.IO.val,2)
                IOtypes(:,i)=cellfun(@(x)sprintf('%s-%d',x,i),IOtypes(:,i),...
                                     'uniformoutput',false);
            end
        else
            % Mixed-variant.
            for j=1:size(s.IO.val,2)
                for i=1:size(s.IO.val,1)
                    IOtypes{i,j}=sprintf('%s-%d(%d)',IOtypes{i,j},j,...
                                         s.IO.struct.block(i,j));
                end
            end
        end
    end
end

if strcmp(sel,'EO') || strcmp(sel,'struct')
    % Set EO parameter types.
    EOtypes={'EX','EY','EZ','om','ph','ka'}';
    if size(s.EO.val,2)>1
        EOtypes=repmat(EOtypes,1,size(s.EO.val,2));
        % Only specify camera ids if any differ from camera sequence number.
        useCamIds=any(1:length(s.EO.id)~=s.EO.id);
        for i=1:size(s.EO.val,2)
            % Id for this camera.
            if useCamIds
                camStr=sprintf('-%d(%d)',i,s.EO.id(i));
            else
                camStr=sprintf('-%d',i);
            end
            EOtypes(:,i)=cellfun(@(x)[x,camStr],EOtypes(:,i),...
                                 'uniformoutput',false);
        end
    end
end

if strcmp(sel,'OP') || strcmp(sel,'struct')
    OPtypes={'OX','OY','OZ'}';
    if size(s.OP.val,2)>1
        OPtypes=repmat(OPtypes,1,size(s.OP.val,2));
        if any(s.prior.OP.isCtrl)
            OPtypes(:,s.prior.OP.isCtrl)=repmat({'CX','CY','CZ'}',1,nnz(s.prior.OP.isCtrl));
        end
        if any(s.prior.OP.isCheck)
            OPtypes(:,s.prior.OP.isCheck)=repmat({'HX','HY','HZ'}',1,nnz(s.prior.OP.isCheck));
        end
        
        for i=1:size(s.OP.val,2)
            % Id for this OP.
            OPstr=sprintf('-%d',i);
            if s.OP.id(i)~=i
                OPstr=[OPstr,sprintf('/%d',s.OP.id(i))];
            end
            if s.OP.rawId(i)~=s.OP.id(i)
                OPstr=[OPstr,sprintf('/%d',s.OP.rawId(i))];
            end
            if ~isempty(s.OP.label{i})
                OPstr=[OPstr,'-',s.OP.label{i}];
            end
            OPtypes(:,i)=cellfun(@(x)[x,OPstr],OPtypes(:,i),...
                                 'uniformoutput',false);
        end
    end
end

% Determine what parameter(s) to return.
switch sel
  case 'struct'
    s.IO.type=IOtypes;
    s.EO.type=EOtypes;
    s.OP.type=OPtypes;
    ret=s;
  case 'IO'
    ret=IOtypes;
  case 'EO'
    ret=EOtypes;
  case 'OP'
    ret=OPtypes;
end
