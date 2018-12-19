function ss=writestats(s,fName,desc)
%WRITESTATS Write project statistics file.
%
%   WRITESTATS(S,FILENAME,DESCRIPTION) writes statistics computed from
%   the DBAT structure S to the file FILENAME. The DESCRIPTION string
%   is written at the beginning of the file.
%
%   The ray angles are computed if they are not available in S.
%   S=WRITESTATS(...) returns the updated structure.
%
%   The computed statistics are:
%       - Totals:
%           - Number of cameras.
%           - Number of images.
%           - Number of control points.
%           - Number of object points.
%           - Number of image points.
%       - Cameras:
%           - Ray count, min, max, mean, histogram, list of lowest.
%           - Angle min, max, mean, histogram, list of lowest.
%       - Control points:
%           - Ray count, min, max, mean, histogram, list of lowest.
%           - Angle min, max, mean, histogram, list of lowest.
%       - Object points:
%           - Ray count, min, max, mean, histogram, list of lowest.
%           - Angle min, max, mean, histogram, list of lowest.

% Function to compute the number of digits in a number.
Digits=@(x) floor(log10(x))+1;
DfmtR=@(n) sprintf('%%%dd',n);
DfmtL=@(n) sprintf('%%-%dd',n);
FfmtR=@(n) sprintf('%%%d.1f',n+2);
SfmtL=@(n) sprintf('%%-%ds',n);
SfmtR=@(n) sprintf('%%%ds',n);

[fid,message]=fopen(fName,'wt+');

fprintf(fid,'%s\n',desc);

fprintf(fid,'\nProject file: %s\n',s.proj.fileName);

fprintf(fid,'\nExecution time stamp: %04d-%02d-%02d %02d:%02d:%02d\n', ...
        floor(clock));

nCp=nnz(s.prior.OP.isCtrl);
fprintf(fid,'\nTotal # OP          : %d\n',size(s.IP.vis,1)-nCp);
fprintf(fid,'Total # CP          : %d\n',nCp);
fprintf(fid,'Total # cams        : %d\n',size(s.IP.vis,2));
fprintf(fid,'Total # image marks : %d\n',nnz(s.IP.vis));
fprintf(fid,'Project units       : %s\n',s.proj.objUnit);

fprintf(fid,'\nProject images: no (id), shortened label, name:\n');
camNoDigits=Digits(size(s.IP.vis,2));
camIdDigits=Digits(max(s.EO.id));
imLabels=s.EO.label;
labelLen=max(cellfun(@length,imLabels));

% Generate shorter labels if necessary.
if labelLen>8
    % Find longest common prefix.
    
    % Find shortest path.
    shortestPath=min(cellfun(@length,imLabels));
    % Convert to array of char.
    labelMat=char(imLabels);
    % Find first position where strings differ.
    eq=min(labelMat,[],1)==max(labelMat,[],1);
    neqPos=find(~eq,1,'first');
    % How many chars to cut?
    if shortestPath<neqPos
        % Leave 3 chars in shortest path.
        cut=shortestPath-3-1;
    else
        % Leave 3 equal chars...
        cut=neqPos-3-1;
        if cut>0
            %  ...unless one of the last equal chars are non-alphanumeric
            nonAlnum=find(~isstrprop(imLabels{1}(cut+(1:3)),'alphanum'),1,'last');
            if nonAlnum
                cut=cut+nonAlnum;
            end
        end
    end
    if cut>4
        % Remove prefix.
        imLabels=cellfun(@(x)x(cut+1:end),imLabels,'uniformoutput',false);
        % Recompute length of longest prefix.
        labelLen=max(cellfun(@length,imLabels));
    end
end

fmt=sprintf('  %s (%s), %s, %%s\n',DfmtL(camNoDigits),DfmtL(camIdDigits),SfmtL(labelLen));
for i=1:length(s.EO.name)
    fprintf(fid,fmt,i,s.EO.id(i),imLabels{i},fullfile(s.proj.imDir,s.EO.name{i}));
end

fprintf(fid,'\n\nIMAGE STATISTICS\n');

% Image ray count
camRayCount=full(sum(s.IP.vis,1));
camRayDigits=Digits(max(camRayCount));

fprintf(fid,'\nImage ray count:\n');
fmt=DfmtR(camRayDigits);
fprintf(fid,['  min : ',fmt,'\n'],min(camRayCount));
fprintf(fid,['  max : ',fmt,'\n'],max(camRayCount));
fprintf(fid,['  mean: ',fmt,'\n'],round(mean(camRayCount)));

[n,edges]=histcounts(camRayCount,'binmethod','sturges');
edges(end)=edges(end)+1;
countDigits=Digits(max(edges));

fprintf(fid,'\nImage with lowest ray count: cam no (id), label, count\n');
[count,i]=sort(camRayCount);
fmt=sprintf('  %s (%s), %s, %s\n',DfmtR(camNoDigits),DfmtR(camIdDigits),SfmtL(labelLen),DfmtR(countDigits));
for j=1:nnz(count<count(min(3,end))*1.1+0.1)
    fprintf(fid,fmt,i(j),s.EO.id(i(j)),imLabels{i(j)},count(j));
end

fprintf(fid,'\nImage ray count histogram: nRays, nCams\n');
fmt=sprintf('  %s-%s: %%d\n',DfmtR(countDigits),DfmtR(countDigits));
for i=1:length(n)
    fprintf(fid,fmt,edges(i),edges(i+1)-1,n(i));
end


% Image angles
if isfield(s,'camRayAng') && ~isempty(s.camRayAng)
    camRayAng=s.camRayAng;
else
    camRayAng=camangles(s,'Computing image ray angles')*180/pi;
    s.camRayAng=camRayAng;
end

fprintf(fid,'\nImage ray angles (deg):\n');
fprintf(fid,'  min : %4.1f\n',min(camRayAng));
fprintf(fid,'  max : %4.1f\n',max(camRayAng));
fprintf(fid,'  mean: %4.1f\n',mean(camRayAng));

fprintf(fid,'\nSmallest image ray angles: cam no (id), label, nRays, angle\n');
fmt=sprintf('  %s (%s), %s, %s, %%4.1f\n',DfmtR(camNoDigits),...
            DfmtR(camIdDigits),SfmtL(labelLen),DfmtR(countDigits));
[ang,i]=sort(camRayAng);
for j=1:nnz(ang<ang(min(3,end))*1.1+0.1)
    fprintf(fid,fmt,i(j),s.EO.id(i(j)),imLabels{i(j)},camRayCount(i(j)),...
            ang(j));
end

fprintf(fid,'\nImage ray angle histogram: angle, count\n');
aa=0:5:90;
aHist=hist(camRayAng,aa);
cDigits=Digits(max(aHist));
fmt=sprintf('  %%2d, %s\n',DfmtR(cDigits));
fprintf(fid,fmt,[aa(:),aHist(:)]');

% Compute ray count for CP + OP.
nRays=full(sum(s.IP.vis,2));
CPix=find(s.prior.OP.isCtrl);
OPix=find(~s.prior.OP.isCtrl);

% CP/OP ray angles
if isfield(s,'rayAng') && ~isempty(s.rayAng)
    rayAng=s.rayAng;
else
    rayAng=angles(s,'Computing CP/OP ray angles')*180/pi;
    s.rayAng=rayAng;
end

ixx={CPix,OPix};
titles={'CONTROL POINT STATISTICS','OBJECT POINT STATISTICS'};
strs={'CP','OP'};

for ii=1:length(ixx)
    ix=ixx{ii};
    
    fprintf(fid,'\n\n%s\n',titles{ii});
    
    if ~any(ix)
        continue;
    end
    
    rayHist=ihist(nRays(ix)+1);

    fprintf(fid,'\n%s ray count:\n',strs{ii});

    rayD=Digits(max(nRays(ix)));

    fprintf(fid,['  min : ',DfmtR(rayD),'\n'],min(nRays(ix)));
    fprintf(fid,['  max : ',DfmtR(rayD),'\n'],max(nRays(ix)));
    fprintf(fid,['  mean: ',FfmtR(rayD),'\n'],mean(nRays(ix)));

    fprintf(fid,'\n%s ray count histogram: nRays, count\n',strs{ii});
    [n,~,count]=find(rayHist);
    fmt=sprintf('  %s, %s\n',DfmtR(Digits(max(n))),DfmtR(Digits(max(count))));
    for i=1:length(n)
        fprintf(fid,fmt,n(i)-1,count(i));
    end
    
    [rays,i]=sort(nRays(ix));
    cut=nnz(rays<rays(min(3,end))*1.1+0.1);
    labelLen=max(cellfun(@length,s.OP.label(ix(i(1:cut)))));
    if labelLen==0
        sTitle='';
        fmt=sprintf('  %s (%s), %s, (%%s)\n',...
                    DfmtR(Digits(max(ix(i(1:cut))))),...
                    DfmtR(Digits(max(s.OP.rawId(ix(i(1:cut)))))),...
                    DfmtR(Digits(rays(cut))));
    else
        sTitle='label, ';
        fmt=sprintf('  %s (%s), %s, %s, (%%s)\n',...
                    DfmtR(Digits(max(ix(i(1:cut))))),...
                    DfmtR(Digits(max(s.OP.rawId(ix(i(1:cut)))))),...
                    SfmtL(labelLen),DfmtR(Digits(rays(cut))));
    end
    fprintf(fid,['\n%s with lowest ray count: %s no (id), %snRays, ' ...
                 '(images with rays)\n'],strs{ii},strs{ii},sTitle);
    for j=1:cut
        imNo=find(s.IP.vis(ix(i(j)),:));
        imStr=strjoin(imLabels(imNo),', ');
        if labelLen==0
            fprintf(fid,fmt,ix(i(j)),s.OP.rawId(ix(i(j))),rays(j),imStr);
        else
            fprintf(fid,fmt,ix(i(j)),s.OP.rawId(ix(i(j))),...
                    s.OP.label{ix(i(j))},rays(j),imStr);
        end
    end
    
    fprintf(fid,'\n%s ray angles:\n',strs{ii});
    fprintf(fid,'  min : %4.1f\n',min(rayAng(ix)));
    fprintf(fid,'  max : %4.1f\n',max(rayAng(ix)));
    fprintf(fid,'  mean: %4.1f\n',mean(rayAng(ix)));

    [ang,i]=sort(rayAng(ix));
    cut=nnz(ang<ang(min(3,end))*1.1+0.1);
    labelLen=max(cellfun(@length,s.OP.label(ix(i(1:cut)))));
    if labelLen==0
        sTitle='';
        fmt=sprintf('  %s (%s), %s, %%4.1f, (%%s)\n',...
                    DfmtR(Digits(max(ix(i(1:cut))))),...
                    DfmtR(Digits(max(s.OP.rawId(ix(i(1:cut)))))),...
                    DfmtR(Digits(max(nRays(ix(i(1:cut)))))));
    else
        sTitle='label, ';
        fmt=sprintf('  %s (%s), %s, %s, %%4.1f, (%%s)\n',...
                    DfmtR(Digits(max(ix(i(1:cut))))),...
                    DfmtR(Digits(max(s.OP.rawId(ix(i(1:cut)))))),...
                    SfmtL(labelLen),DfmtR(Digits(max(nRays(ix(i(1:cut)))))));
    end
    fprintf(fid,['\nSmallest %s ray angles: %s no (id), %snRays, ' ...
                 'angle, (images with rays)\n'],strs{ii},strs{ii},sTitle);
    for j=1:cut
        imNo=find(s.IP.vis(ix(i(j)),:));
        imStr=strjoin(imLabels(imNo),', ');
        if labelLen==0
            fprintf(fid,fmt,ix(i(j)),s.OP.rawId(ix(i(j))),...
                    nRays(ix(i(j))),ang(j),imStr);
        else
            fprintf(fid,fmt,ix(i(j)),s.OP.rawId(ix(i(j))),...
                    s.OP.label{ix(i(j))},nRays(ix(i(j))),...
                    ang(j),imStr);
        end
    end
    
    fprintf(fid,'\n%s ray angle histogram: angle, count\n',strs{ii});
    aa=0:5:90;
    aHist=hist(rayAng(ix),aa);
    cDigits=Digits(max(aHist));
    fmt=sprintf('  %%2d, %s\n',DfmtR(cDigits));
    fprintf(fid,fmt,[aa(:),aHist(:)]');
end

fclose(fid);

if nargout>0
    ss=s;
end
