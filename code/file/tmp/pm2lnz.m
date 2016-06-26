OP=round(prob.objPts(:,2:3)'*7);
OPid=prob.objPts(:,1)';

for im=1:length(prob.images)
    fprintf('\n\nImage %d: %s\n',im,prob.images(im).imName);


    ptsXY=prob.markPts(prob.markPts(:,1)==im-1,3:4)';
    ptsId=prob.markPts(prob.markPts(:,1)==im-1,2)';

    for i=1:length(ptsId)
        xy=ptsXY(:,i);
        j=find(ptsId(i)==OPid);
        if j
            fprintf(['      <corner img_x="%.2f" img_y="%.2f" obj_x="%.2f" ' ...
                     'obj_y="%.2f" valid="true"/>\n'],xy,OP(:,j));
        end
    end
end
