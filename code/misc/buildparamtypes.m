function ret=buildparamtypes(s,sel)
%BUILDPARAMTYPES Build parameter type strings.
%
%   S1=BUILDPARAMTYPES(S0) populates the paramTypes fields of the DBAT
%   struct S1 based on the information in the DBAT struct S0, as
%   described below.
%
%   T=BUILDPARAMTYPES(S,SEL), where SEL is 'IO', 'EO', or 'OP',
%   returns the cell array of selected parameter types instead.
%
%   Each S.IO column stores the parameters below. The first 10 may
%   be estimated by the bundle.
%       px,
%       py      - principal point in camera units (typically mm).
%       cc      - camera constant in camera units.
%       K1,
%       K2,
%       K3      - radial distortion parameters of Brown (1971).
%       P1,
%       P2      - tangential distortion parameters of Brown (1971).
%       fa      - affine parameter. Aspect will be (1+fa):1.
%       fs      - skew parameters.
%       sw,
%       sh      - sensor width and height in camera units.
%       iw,
%       ih      - image width in pixels.
%       rx,
%       ry      - image resolution.
%
%   The names above are returned in the 16-by-nCams IO cell array. If
%   multiple S.IO columns are present, the column number is
%   appended. If the block number is different from the image
%   number, the block number is appended in parantheses.
%
%   Each S.EO column stores the parameters below. The first 6 parameters
%   may be estimated by the bundle.
%       EX,
%       EY,
%       EZ       - external coordinates of the camera center in project units.
%       omega,
%       phi,
%       kappa   - Euler angles for the camera orientation in radians.
%       tt      - parameter indicating which Euler convention is
%                 used. Currently only t=0 (omega, phi,kappa) is supported.
%
%   The first two letters of the names above are returned in the
%   7-by-nImages EO cell array. If multiple EO columns are present, a
%   camera identifier is appended to each parameter. The camera
%   parameter is consists of the camera sequence number and camera id.
%
%   Each S.OP column stores the X, Y, Z coordinates. The OP names are
%   returned in the 3-by-nOP cell array OP. Object points are prefixed
%   with 'O', i.e. 'OX', 'OY', 'OZ'. Control points are prefixed with
%   'C'. Check points are prefixed with 'H'. Furthermore, a point
%   identifier is appended. The point identifier consists of the
%   sequence number and if necessary, the OP id, the OP raw id, and
%   the OP label.
%
%See also: PROB2DBATSTRUCT.

if nargin<2, sel='struct'; end

IOtypes={};
EOtypes={};
OPtypes={};

if strcmp(sel,'IO') || strcmp(sel,'struct')
    % IO parameter types.
    Knames=arrayfun(@(x)sprintf('K%d',x),1:s.nK,'uniformoutput',false);
    Pnames=arrayfun(@(x)sprintf('P%d',x),1:s.nP,'uniformoutput',false);
    IOtypes={'px','py','cc',Knames{:},Pnames{:},'fa','fs','sw', ...
             'sh','iw','ih','rx','ry'}';
    if size(s.IO,2)>1
        IOtypes=repmat(IOtypes,1,size(s.IO,2));
        if any(any(s.IOblock~=repmat(1:size(s.IO,2),size(s.IO,1),1)))
            % At least partly block-invariant.
            for j=1:size(s.IO,2)
                for i=1:size(s.IO,1)
                    IOtypes{i,j}=sprintf('%s-%d(%d)',IOtypes{i,j},j,...
                                         s.IOblock(i,j));
                end
            end
        else
            % Purely image-invariant.
            for i=1:size(s.IO,2)
                IOtypes(:,i)=cellfun(@(x)sprintf('%s-%d',x,i),IOtypes(:,i),...
                                     'uniformoutput',false);
            end
        end
    end
end

if strcmp(sel,'EO') || strcmp(sel,'struct')
    % Set EO parameter types.
    EOtypes={'EX','EY','EZ','om','ph','ka','tt'}';
    if size(s.EO,2)>1
        EOtypes=repmat(EOtypes,1,size(s.EO,2));
        % Only specify camera ids if any differ from camera sequence number.
        useCamIds=any(1:length(s.camIds)~=s.camIds);
        for i=1:size(s.EO,2)
            % Id for this camera.
            if useCamIds
                camStr=sprintf('-%d(%d)',i,s.camIds(i));
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
    if size(s.OP,2)>1
        OPtypes=repmat(OPtypes,1,size(s.OP,2));
        if any(s.isCtrl)
            OPtypes(:,s.isCtrl)=repmat({'CX','CY','CZ'}',1,nnz(s.isCtrl));
        end
        if any(s.isCheck)
            OPtypes(:,s.isCheck)=repmat({'HX','HY','HZ'}',1,nnz(s.isCheck));
        end
        
        for i=1:size(s.OP,2)
            % Id for this OP.
            OPstr=sprintf('-%d',i);
            if s.OPid(i)~=i
                OPstr=[OPstr,sprintf('/%d',s.OPid(i))];
            end
            if s.OPrawId(i)~=s.OPid(i)
                OPstr=[OPstr,sprintf('/%d',s.OPrawId(i))];
            end
            if ~isempty(s.OPlabels{i})
                OPstr=[OPstr,'-',s.OPlabels{i}];
            end
            OPtypes(:,i)=cellfun(@(x)[x,OPstr],OPtypes(:,i),...
                                 'uniformoutput',false);
        end
    end
end

% Determine what parameter(s) to return.
switch sel
  case 'struct'
    s.paramTypes=struct('IO',{IOtypes},'EO',{EOtypes},'OP',{OPtypes});
    ret=s;
  case 'IO'
    ret=IOtypes;
  case 'EO'
    ret=EOtypes;
  case 'OP'
    ret=OPtypes;
end
