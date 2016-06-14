% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));
problem=2

dampings={'none','gna','lm','lmp'};

dampings=dampings(2);

if ~exist('fName','var')
    switch problem
      case 1
        stub='fixed1-5';
      case 2
        stub='weighted1-5';
      case 3
        stub='w15-op1';
      case 4
        stub='wsmart-no-orient';
      otherwise
        error('bad problem number');
    end
    dataDir=stub;
    fName=fullfile(curDir,'data','weighted','pm','sxb',dataDir,...
                   [stub,'-pmexport.txt']);
    cpName=fullfile(curDir,'data','weighted','pm','sxb',dataDir,'ctrlpts.txt');
    ptName=fullfile(curDir,'data','weighted','pm','sxb',dataDir,[stub,'-3dpts.txt']);
    fprintf('No data file specified, using ''%s''.\n',fName);
    disp(['Set variable ''fName'' to name of Photomodeler Export file if ' ...
          'you wish to use another file.']);
    disp(' ')
end

if ~exist('prob','var')
    fprintf('Loading data file %s...',fName);
    prob=loadpm(fName);
    prob0=prob;
    fprintf('done.\n');
    if any(isnan(cat(2,prob.images.imSz)))
        error('Image sizes unknown!');
    end
    disp('done.')
    ctrlPts=loadcpt(cpName);
    fprintf('done.\n');
    if true
        if ~all(ismember(prob.ctrlPts(:,1),ctrlPts.id))
            error('Control point id mismatch.');
        end
        fprintf('Loading 3D point table %s...',ptName);
        pts=loadpm3dtbl(ptName);
        fprintf('done.\n');
    
        % Determine offset between real positions and positions used in
        % the bundle.
        [~,ia,ib]=intersect(prob.ctrlPts(:,1),pts.id);
        offset=prob.ctrlPts(ia,2:4)'-pts.pos(:,ib);

        % Offset range should be as small as difference between the number of
        % digits used. Warn if it is larger than 1e-3 object units
        % (typically 1mm).
        offsetRange=max(offset,[],2)-min(offset,[],2);
        if max(offsetRange)>1e-3
            warning('Large offset range:')
            offsetRange
        end

        % Adjust a priori control point positions by the offset.
        ctrlPts.pos=ctrlPts.pos+repmat(offset,1,size(ctrlPts,2));
    
        % Replace a posteriori ctrl positions and std by a priori values.
        [~,ia,ib]=intersect(ctrlPts.id,prob.ctrlPts(:,1));
        prob.ctrlPts(ib,2:4)=ctrlPts.pos(:,ia)';
        prob.ctrlPts(ib,5:7)=ctrlPts.std(:,ia)';
    
        cpId=prob.ctrlPts(:,1);
    end
else
    disp('Using pre-loaded data. Do ''clear prob'' to reload.');
end
ss0=prob2dbatstruct(prob);

% Optionally, change estimate/prior handling.

s0=ss0;

% Warn for non-uniform mark std.
uniqueSigmas=unique(s0.markStd(:));

if length(uniqueSigmas)~=1
    warning('Multiple sigmas, assuming sigma=1');
    s0.prior.sigmas(1)=1;
else
    s0.prior.sigmas=uniqueSigmas;
end

% Reset all parameters to be estimated.
s0.IO(s0.estIO)=nan;
s0.IO(s0.useIOobs)=s0.prior.IO(s0.useIOobs);

s0.EO(s0.estEO)=nan;
s0.EO(s0.useEOobs)=s0.prior.EO(s0.useEOobs);

s0.OP(s0.estOP)=nan;
s0.OP(s0.useOPobs)=s0.prior.OP(s0.useOPobs);

% Use specified sigma as first approximation.
s0.markStd(:)=s0.prior.sigmas(1);

s1=resect(s0,'all',cpId,1,0,cpId);
s2=forwintersect(s1,'all',true);

s2.x0desc='Camera calibration from EXIF value';

s=s2;
h=plotnetwork(s,'title','Initial network',...
              'axes',tagfigure(sprintf('network%d',i)),...
              'camsize',0.1);

if printdemofigures
    h=get(h,'parent');
    figDir=fullfile('..','docsrc','manual','ill');
    files={'ccamx0.eps'};
    for i=1:length(h)
        print(h(i),'-depsc2',fullfile(figDir,files{i}));
    end
end

result=cell(size(dampings));
ok=nan(size(dampings));
iters=nan(size(dampings));
sigma0=nan(size(dampings));
E=cell(size(dampings));

for i=1:length(dampings)
    fprintf('Running the bundle with damping %s...\n',dampings{i});

    % Run the bundle.
    [result{i},ok(i),iters(i),sigma0(i),E{i}]=bundle(s,dampings{i},'trace','dofverb');
    
    if ok(i)
        fprintf('Bundle ok after %d iterations with sigma0=%.2f (%.2f pixels)\n', ...
                iters(i),sigma0(i),sigma0(i)*s.prior.sigmas(1));
    else
        fprintf(['Bundle failed after %d iterations. Last sigma0 estimate=%.2f ' ...
                 '(%.2f pixels)\n'],iters(i),sigma0(i),...
                sigma0(i)*s.prior.sigmas(1));
    end
end

resFile=strrep(fName,'.txt','_result_file.txt');

COP=bundle_result_file(result{1},E{1},resFile);
OPstd=full(reshape(sqrt(diag(COP)),3,[]));
CEO=bundle_cov(result{1},E{1},'CEO');
EOstd=reshape(full(sqrt(diag(CEO))),6,[]);
EOposStd=EOstd(1:3,:);
EOangStd=EOstd(4:6,:)*180/pi;

fprintf('\nBundle result file %s generated.\n',resFile);

h=plotparams(result{1},E{1});

if printdemofigures
    figDir=fullfile('..','docsrc','manual','ill');
    files={'ccamiotrace.eps','ccameotrace.eps','ccamoptrace.eps', ...
           'ccamgnatrace.eps'};
    for i=1:length(h)
        print(h(i),'-depsc2',fullfile(figDir,files{i}));
    end
end

h=plotcoverage(result{1},true);

if printdemofigures
    figDir=fullfile('..','docsrc','manual','ill');
    files={'ccamcoverage.eps'};
    for i=1:length(h)
        print(h(i),'-depsc2',fullfile(figDir,files{i}));
    end
end

h=plotimagestats(result{1},E{1});

if printdemofigures
    figDir=fullfile('..','docsrc','manual','ill');
    files={'ccamimstats.eps'};
    for i=1:length(h)
        print(h(i),'-depsc2',fullfile(figDir,files{i}));
    end
end

h=plotopstats(result{1},E{1},COP);

if printdemofigures
    figDir=fullfile('..','docsrc','manual','ill');
    files={'ccamopstats.eps'};
    for i=1:length(h)
        print(h(i),'-depsc2',fullfile(figDir,files{i}));
    end
end

if printdemofigures, doPause=0; else doPause='on'; end

for i=1:length(E)
    h=plotparams(result{i},E{i},'noio','noeo','noop');
    fig=tagfigure(sprintf('network%d',i));
    fprintf('Displaying bundle iteration playback for method %s in figure %d.\n',E{i}.damping.name,double(fig));
    h=plotnetwork(result{i},E{i},'title',...
                  ['Damping: ',dampings{i},'. Iteration %d of %d'], ...
                  'axes',fig,'pause',doPause,'camsize',0.1); 
end

if printdemofigures
    h=get(h,'parent');
    figDir=fullfile('..','docsrc','manual','ill');
    files={'ccamxfinal.eps'};
    for i=1:length(h)
        print(h(i),'-depsc2',fullfile(figDir,files{i}));
    end
end
