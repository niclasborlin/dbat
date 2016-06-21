function romabundledemo
%ROMABUNDLEDEMO Bundle demo for DBAT.
%
%   ROMABUNDLEDEMO runs the bundle on the PhotoModeler export file
%   of the ROMA data set. The PhotoModeler EO values are used as
%   initial values, except that the EO position are disturbed by
%   random noise with sigma=0.1 m. The OP initial values are
%   computed by forward intersection.

% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

% Base dir with input files.
inputDir=fullfile(curDir,'data','phor2013');

% PhotoModeler text export file and report file.
inputFile=fullfile(inputDir,'pmexports','roma-pmexport.txt');
% Report file name.
reportFile=fullfile(inputDir,'dbatexports','roma-dbatreport.txt');;

fprintf('Loading data file %s...',inputFile);
prob=loadpm(inputFile);
probRaw=prob;
if any(isnan(cat(2,prob.images.imSz)))
    error('Image sizes unknown!');
end
disp('done.')

% Convert loaded PhotoModeler data to DBAT struct.
s0=prob2dbatstruct(prob);
ss0=s0;

% Don't estimate IO data, i.e. treat it as exact.  This block is not
% really necessary, but may serve as a starting point if IO
% parameters are to be estimated.
s0.IO=s0.prior.IO;
s0.estIO=false(size(s0.IO));
s0.useIOobs=false(size(s0.IO));

% Noise sigma [m].
noiseLevel=0.1;

% Use supplied EO data as initial values. Again, this block is not
% really necessary but may serve as a starting point for modifications.
s0.EO=s0.prior.EO;
s0.estEO(1:6,:)=true; % 7th element is just the axis convention.
s0.useEOobs=false(size(s0.EO));
s0.EO(1:3,:)=s0.EO(1:3,:)+randn(3,size(s0.EO,2))*noiseLevel;

% Copy CP values and treat them as fixed.
s0.OP(:,s0.isCtrl)=s0.prior.OP(:,s0.isCtrl);
s0.estOP=repmat(~s0.isCtrl(:)',3,1);
s0.useOPobs=repmat(s0.isCtrl(:)',3,1);
% Compute initial OP values by forward intersection.
correctedPt=reshape(pm_multilenscorr1(diag([1,-1])*s0.markPts,s0.IO,3,2,...
                                      s0.ptCams,size(s0.IO,2)),2,[]);
s0.OP(:,~s0.isCtrl)=pm_multiforwintersect(s0.IO,s0.EO,s0.cams,s0.colPos,correctedPt,find(~s0.isCtrl));

% Warn for non-uniform mark std.
uniqueSigmas=unique(s0.markStd(:));

if length(uniqueSigmas)~=1
    uniqueSigmas
    error('Multiple mark point sigmas')
end

if all(uniqueSigmas==0)
    s0.markStd(:)=1;
end

% Fix the datum by fixing camera 1...
s0.estEO(:,1)=false;
% ...and the largest other absolute camera coordinate.
camDiff=abs(s0.EO(1:3,:)-repmat(s0.EO(1:3,1),1,size(s0.EO,2)));
[i,j]=find(camDiff==max(camDiff(:)));
s0.estEO(i,j)=false;

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
        fprintf('Bundle ok after %d iterations with sigma0=%.2f pixels\n', ...
                iters(i),sigma0(i));
    else
        fprintf(['Bundle failed after %d iterations. Last sigma0 estimate=%.2f ' ...
                 'pixels\n'],iters(i),sigma0(i));
    end
end

COP=bundle_result_file(result{1},E{1},reportFile);

fprintf('\nBundle result file %s generated.\n',resFile);

plotparams(result{1},E{1},'noop');

plotcoverage(result{1},true);

plotimagestats(result{1},E{1});

plotopstats(result{1},E{1},COP);

for i=1:length(E)
    fig=tagfigure(sprintf('network%d',i));
    fprintf('Displaying bundle iteration playback for method %s in figure %d.\n',E{i}.damping.name,double(fig));
    plotnetwork(result{i},E{i},'trans','up','align',1,'title',...
                ['Damping: ',dampings{i},'. Iteration %d of %d'], ...
                'axes',fig,'pause','on');
end
