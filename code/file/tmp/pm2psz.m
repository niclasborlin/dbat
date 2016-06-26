% Control points.

cpIds=1001:1004;

for i=1:4
    id=1000+i;
    j=find(prob.objPts(:,1)==id);
    if any(j)
        fprintf('        <marker id="%d" label="CP%d">\n',i,i);
        fprintf('          <reference x="%.16e" y="%.16e" z="%.16e" enabled="true"/>\n',...
                prob.objPts(j,2:4));
        fprintf('        </marker>\n');
    end
end

fprintf('\n\n');

fprintf('          <markers>\n');

for i=1:4
    id=1000+i;
    if any(prob.markPts(:,2)==id)
        fprintf('            <marker marker_id="%d">\n',i);
        ims=sort(prob.markPts(prob.markPts(:,2)==id,1));
        for k=1:length(ims)
            xy=prob.markPts(prob.markPts(:,1)==ims(k) & ...
                            prob.markPts(:,2)==id,3:4);
            fprintf('              <location camera_id="%d" pinned="true" x="%.7e" y="%.7e"/>\n',ims(k),xy);
        end
        fprintf('            </marker>\n');
    end
end

fprintf('          </markers>\n');

tieIds=setdiff(intersect(prob.markPts(:,2),prob.objPts(:,1)),cpIds);
newIds=sparse(tieIds,1,0:length(tieIds)-1);

for i=1:length(prob.images)
    pts=prob.markPts(prob.markPts(:,1)==i-1,:);
    pts=pts(ismember(pts(:,2),tieIds),:);
    [~,j]=sort(pts(:,2));
    pts=pts(j,2:4);
    clear Data
    Data.vertex.id=full(newIds(pts(:,1)));
    Data.vertex.x=pts(:,2);
    Data.vertex.y=pts(:,3);
    ply_write(Data,sprintf('/tmp/projections%d.ply',i-1),'binary_little_endian');
end
