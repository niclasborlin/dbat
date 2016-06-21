% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

stub='weighted-shifted';
fName=fullfile(curDir,'data','weighted','ps','sxb',stub,[stub,'.psz']);
    
s=loadpsz2(fName);

imSz=s.camera.imSz(:);
defCam=[s.camera.focal;s.camera.pp(:);s.camera.sensorFormat(:);zeros(5,1)];

job=struct('title','Photoscan import','defCam',defCam,'defCamStd',zeros(size(defCam)),'imSz',imSz);

% Camera positions.
CC=s.global.CC;
% Angles
ang=nan(size(CC));
for i=1:size(ang,2)
    RR=s.global.R(:,:,i);
    ang(:,i)=derotmat3d(RR)';
end
angPM=ang([3,2,1],:)*180/pi;

images=repmat(struct('imName','','outer',nan(1,6),'outerStd',zeros(1,6),...
                     'imSz',imSz),size(CC,2),1);
for i=1:length(images)
    images(i).imName=s.imNames{i};
    images(i).outer=[CC(:,i);angPM(:,i)]';
end

ctrlPts=s.global.ctrlPts;

objPts=[ctrlPts;[s.global.objPts,nan(size(s.global.objPts,1),3)]];

markPts=[s.markPts.all,1*ones(size(s.markPts.all,1),2)];

prob2=struct('job',job,'images',images,'ctrlPts',ctrlPts,'objPts',objPts,...
             'markPts',markPts);

s0=prob2dbatstruct(prob2);

% Don't estimate IO data (treat it as exact).
s0.IO=s0.prior.IO; % Not really necessary...
s0.estIO=false(size(s0.IO));
s0.useIOobs=false(size(s0.IO));
%s0.estIO(1:3)=true;

% Use supplied EO data as initial values. Treat EO data as free.
s0.EO=s0.prior.EO; % No really necessary...
s0.estEO(1:6,:)=true;
s0.useEOobs=false(size(s0.EO));

% Use supplied OP data as initial values. Treat control points as
% exact.
s0.OP(:,s0.isCtrl)=s0.prior.OP(:,s0.isCtrl);
s0.estOP=true(size(s0.OP));
s0.useOPobs=repmat(s0.isCtrl',3,1);

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
    [result{i},ok(i),iters(i),sigma0(i),E{i}]=bundle(s0,dampings{i},'trace');
    
    if ok(i)
        fprintf('Bundle ok after %d iterations with sigma0=%.2f (%.2f pixels)\n', ...
                iters(i),sigma0(i),sigma0(i)*s0.prior.sigmas(1));
    else
        fprintf(['Bundle failed after %d iterations. Last sigma0 estimate=%.2f ' ...
                 '(%.2f pixels)\n'],iters(i),sigma0(i),...
                sigma0(i)*s0.prior.sigmas(1));
    end
end

[p,n,e]=fileparts(fName);

resFile=fullfile(p,[n,'_result_file.txt']);

COP=bundle_result_file(result{1},E{1},resFile);
OPstd=full(reshape(sqrt(diag(COP)),3,[]));
CEO=bundle_cov(result{1},E{1},'CEO');
EOstd=reshape(full(sqrt(diag(CEO))),6,[]);
EOposStd=EOstd(1:3,:);
EOangStd=EOstd(4:6,:)*180/pi;

fprintf('\nBundle result file %s generated.\n',resFile);

plotparams(result{1},E{1},'noop');

plotcoverage(result{1},true);

plotimagestats(result{1},E{1});

plotopstats(result{1},E{1},COP);

for i=1:length(E)
    fig=tagfigure(sprintf('network%d',i));
    fprintf('Displaying bundle iteration playback for method %s in figure %d.\n',E{i}.damping.name,double(fig));
    plotnetwork(result{i},E{i},'title',...
                ['Damping: ',dampings{i},'. Iteration %d of %d'], ...
                'axes',fig,'camsize',100,'pause','on');
end
