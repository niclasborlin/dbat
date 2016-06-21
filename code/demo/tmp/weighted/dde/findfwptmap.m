% Load tables
tbl90=load(sprintf('edfdata/orig/%d_cpids.txt',1990));
tbl07=load(sprintf('edfdata/orig/%d_cpids.txt',2007));

% Find out which levels there are
lvl90=unique(-tbl90(tbl90(:,2)<0,2));
lvl07=unique(-tbl07(tbl07(:,2)<0,2));
lvls=intersect(lvl07,lvl90);

% For each level...
for i=1:length(lvls)
    % Extract points from each year
    pt90=tbl90(tbl90(:,2)==-lvls(i),[1,4:6]);
    pt07=tbl07(tbl07(:,2)==-lvls(i),[1,4:6]);
    % Calculate distances between all points.
    D=zeros(size(pt90,1),size(pt07,1));
    for j=1:size(pt90,1)
        for k=1:size(pt07,1)
            D(j,k)=norm(pt90(j,2:4)-pt07(k,2:4));
        end
    end
    id90=pt90(:,1);
    id07=pt07(:,1);
    if size(D,1)~=size(D,2)
        D(max(size(D)),max(size(D)))=0;
        if length(id90)<size(D,1)
            id90(size(D,1))=0;
        end
        if length(id07)<size(D,2)
            id07(size(D,2))=0;
        end
    end
    % Find assignment that gives total cost (distance).
    C=hungarian(D);
    % Remap cost matrix and row id.
    D=D(C,:);
    id90=id90(C);
    % Extract each distance.
    d=diag(D);
    % OK mappings have a motion <1m
    ok=id90~=0 & id07~=0 & d<1;
    % Non-OK mappings have a motion >1m
    nonOK=id90~=0 & id07~=0 & d>=1;
    disp(sprintf('Found %d good mappings (max=%f) on level %d',...
                 nnz(ok),max(d(ok)),lvls(i)))
    s=sprintf('Found %d bad mappings ',nnz(nonOK));
    if any(nonOK)
        ids=[id90(nonOK),id07(nonOK)];
        s=[s,sprintf('%d-%d, ',ids')];
        s=[s,sprintf('(min=%f) ',min(d(nonOK)))];
    end
    s=[s,sprintf('on level %d',lvls(i))];
    disp(s)
end
