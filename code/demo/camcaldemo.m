% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

% Defult to Olympus Camedia C4040Z dataset if no data file is specified.
if ~exist('fName','var')
    fName=fullfile(curDir,'data','C4040Z-2272x1704.txt');
    fprintf('No data file specified, using ''%s''.\n',fName);
    disp(['Set variable ''fName'' to name of Photomodeler Export file if ' ...
          'you wish to use another file.']);
    disp(' ')
end

if ~exist('prob','var')
    fprintf('Loading data file %s...',fName);
    prob=loadpm(fName);
    disp('done.')
else
    disp('Using pre-loaded data. Do ''clear prob'' to reload.');
end
s0=prob2dbatstruct(prob);

% Fix the datum by fixing CP 1, 2 + the Z coordinate of OP 3.
s0.cOP(:,ismember(s0.OPid,1001:1002))=false;
s0.cOP(3,ismember(s0.OPid,1003))=false;

% Estimate px,py,c,K1-K3,P1-P2.
s0.cIO(1:8,:)=true;

dampings={'none','gna','lm','lmp'};

dampings=dampings(2);

result=cell(size(dampings));
ok=nan(size(dampings));
iters=nan(size(dampings));
sigma0=nan(size(dampings));
E=cell(size(dampings));

for i=1:length(dampings)
    fprintf('Running the bundle with damping %s...\n',dampings{i});

    % Run the bundle.
    [result{i},ok(i),iters(i),sigma0(i),E{i},CXX,CIOF,CIO,CEOF,CEO,COPF,COP]=...
        bundle(s0,dampings{i},'trace','cxx','ciof','cio','ceof','ceo', ...
               'copf','cop');
    
    if true
        disp('Comparing CIOF with CIO')
        disp(full(max(max(abs(mkblkdiag(CIOF,size(s0.cIO,1))-CIO)))))
        disp('Comparing CEOF with CEO')
        disp(full(max(max(abs(mkblkdiag(CEOF,size(s0.cEO,1))-CEO)))))
        disp('Comparing COPF with COP')
        disp(full(max(max(abs(mkblkdiag(COPF,size(s0.cOP,1))-COP)))))

        nIO=nnz(s0.cIO);
        nEO=nnz(s0.cEO);
        nOP=nnz(s0.cOP);
        disp('Comparing CIOF with CXX')
        disp(full(max(max(abs(CIOF(s0.cIO(:),s0.cIO(:))-...
                              CXX(1:nIO,1:nIO))))))
        disp('Comparing CEOF with CXX')
        disp(full(max(max(abs(CEOF(s0.cEO(:),s0.cEO(:))-...
                              CXX(nIO+(1:nEO),nIO+(1:nEO)))))))
        disp('Comparing COPF with CXX')
        disp(full(max(max(abs(COPF(s0.cOP(:),s0.cOP(:))-...
                              CXX(nIO+nEO+(1:nOP),nIO+nEO+(1:nOP)))))))
    end

    if ok(i)
        fprintf('Bundle ok after %d iterations with sigma0=%.2f pixels\n', ...
                iters(i),sigma0(i));
    else
        fprintf(['Bundle failed after %d iterations. Last sigma0 estimate=%.2f ' ...
                 'pixels\n'],iters(i),sigma0(i));
    end
end

% Rotate to have +Z up.
T0=blkdiag(1,[0,-1;1,0],1);

for i=1:length(E)
    plotnetwork(result{i},E{i},'trans',T0,'align',1,'title',...
                ['Damping: ',dampings{i},'. Iteration %d of %d'], ...
                'axes',figure(i),'pause','on','camerasize',0.25);
end
