function [ok,msg]=writepsmarkertbl(psTbl,fName)
%WRITEPSMARKERTBL Write Photoscan marker table to a text file.
%
%    [OK,MSG]=WRITEPSMARKERTBL(PSTBL,FNAME) writes a Photoscan marker
%    table PSTBL to the text file FNAME, corresponding to the
%    <markers> section of the doc.xml file within a Photoscan
%    archive. If the file name is omitted, the text is written to
%    stdout instead. On failure a return status and message is returned.


msg='';

if nargin<2
    fid=1;
else
    [fid,msg]=fopen(fName,'wt');
    if fid<0
        ok=false;
        return;
    end
end

ft={'false','true'};

fprintf(fid,'          <markers>\n');
ids=unique(psTbl.markerId);
for i=1:length(ids)
    fprintf(fid,'            <marker marker_id="%d">\n',ids(i));
    for j=find(psTbl.markerId==ids(i))'
        fprintf('              <location camera_id="%d" pinned="%s" x="%.10f" y="%.10f"/>\n',...
                psTbl.cameraId(j),ft{psTbl.pinned(j)+1},psTbl.pos(:,j));
    end
    fprintf(fid,'            </marker>\n');
end
fprintf(fid,'          </markers>\n');

if fid>2
    fclose(fid);
end
