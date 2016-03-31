runAsBundle=true;
weighted=true;
cpSigma=1e-2;

% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

dampings={'none','gna','lm','lmp'};

dampings=dampings(2);

% Defult to Olympus Camedia C4040Z dataset if no data file is specified.
if ~exist('fName','var')
    if runAsBundle
        if weighted
            switch cpSigma
              case 1e-6
                wStr='1micron';
              case 1e-3
                wStr='1mm';
              case 1e-2
                wStr='1cm';
              otherwise
                error('Bad weight')
            end
            d=['weighted-bundle-',wStr];
            f=['c4040z-',wStr,'-pmexport.txt'];
        else
            d='fixed-as-bundle';
            f='c4040z-bundle-pmexport.txt';
        end
    else
        d='fixed';
        f='c4040z-pmexport.txt';
    end
    fName=fullfile(curDir,'data','weighted','camcal',d,f);
    cpName=fullfile(curDir,'data','weighted','camcal',d,'ctrlpts.txt');
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
    if runAsBundle
        % Inject prior control point information into what's alrady loaded.
        if ~isequal(sort(cpId(:)),sort(prob.ctrlPts(:,1)))
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
    else
        % Add standard control points for camera calibration project.
        if ~isempty(prob.ctrlPts)
            error('CTRL PTS info exists.');
        end
        prob.ctrlPts=[cpId;CP;CPS]';
    end
    disp('done.')
else
    disp('Using pre-loaded data. Do ''clear prob'' to reload.');
end
ss0=prob2dbatstruct(prob);

s0=ss0;

% Warn for non-uniform mark std.
uniqueSigmas=unique(s0.markStd(:));

if length(uniqueSigmas)~=1
    warning('Multiple sigmas, assuming sigma=.1');
    s0.prior.sigmas(1)=0.1;
else
    s0.prior.sigmas=uniqueSigmas;
end

if runAsBundle
    % Do nothing, estIO is already false.
else
    % Estimate px,py,c,K1-K3,P1-P2.
    s0.estIO([1:5,7:8],:)=true;
    % No prior observations.
    s0.useIOobs(:)=false;

    % Set initial IO parameters.
    s0.IO(1)=s0.IO(11)/2;  % px = center of sensor
    s0.IO(2)=-s0.IO(12)/2; % py = center of sensor (sign is due to camera model)
    s0.IO(3)=7.3;          % c = EXIF value.
    s0.IO(4:8)=0;          % K1-K3, P1-P2 = 0.
end
   
% Clear EO and OP parameters.
s0.EO(s0.estEO)=nan;
s0.OP(s0.estOP)=nan;
% Insert any prior obs to use.
s0.EO(s0.useEOobs)=s0.prior.EO(s0.useEOobs);
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
