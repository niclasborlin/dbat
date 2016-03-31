% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

dampings={'none','gna','lm','lmp'};

dampings=dampings(2);

if ~exist('fName','var')
    dataVersion='fixed';
    fName=fullfile(curDir,'data','weighted','sxb',dataVersion,'gcponly-pmexport.txt');
    cpName=fullfile(curDir,'data','weighted','sxb',dataVersion,'ctrlpts.txt');
    fprintf('No data file specified, using ''%s''.\n',fName);
    disp(['Set variable ''fName'' to name of Photomodeler Export file if ' ...
          'you wish to use another file.']);
    disp(' ')
end

if ~exist('prob','var')
    fprintf('Loading data file %s...',fName);
    prob=loadpm(fName);
    if any(isnan(cat(2,prob.images.imSz)))
        error('Image sizes unknown!');
    end
    [cpId,CP,CPS,cpNames]=loadcpt(cpName);
    if ~all(ismember(prob.ctrlPts(:,1),cpId))
        error('Control point id mismatch.');
    end
    % Replace XYZ positions.
    [~,ia,ib]=intersect(cpId,prob.ctrlPts(:,1));
    % Determine offset.
    offset=mean(CP(:,ia)'-prob.ctrlPts(ib,2:4),1)';
    prob.ctrlPts(ib,2:4)=CP(:,ia)'-repmat(offset',length(ia),1);
    % Replace sigmas.
    prob.ctrlPts(ib,5:7)=CPS(:,ia)';
    % Keep track of control point names.
    names=cell(size(cpNames));
    names(ib)=cpNames(ia);
    cpNames=names;
    disp('done.')
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
    [result{i},ok(i),iters(i),sigma0(i),E{i}]=bundle(s,dampings{i},'trace');
    
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
