runAsBundle=true;
weighted=true;

% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

dampings={'none','gna','lm','lmp'};

dampings=dampings(2);

% Defult to Olympus Camedia C4040Z dataset if no data file is specified.
if ~exist('fName','var')
    if runAsBundle
        if weighted
            stub='weighted-bundle-1mm';
        else
            stub='fixed-as-bundle';
        end
        d=stub;
        f=[stub,'-pmexport.txt'];
    else
        stub='fixed';
        d=stub;
        f=[stub,'-pmexport.txt'];
    end
    fName=fullfile(curDir,'data','weighted','pm','camcal',d,f);
    cpName=fullfile(curDir,'data','weighted','pm','camcal',d,'ctrlpts.txt');
    ptName=strrep(fName,'-pmexport.txt','-3dpts.txt');
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
    if runAsBundle
        % Inject prior control point information into what's alrady loaded.
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
        
        prob.ctrlPts(:,5:end)=prob0.ctrlPts(:,5:end);
    else
        % Add standard control points for camera calibration project.
        fprintf('Loading 3D point table %s...',ptName);
        pts=loadpm3dtbl(ptName);
        fprintf('done.\n');
        offset=zeros(3,4);
        if ~isempty(prob.ctrlPts)
            error('CTRL PTS info exists.');
        end
        prob.ctrlPts=[ctrlPts.id;ctrlPts.pos;ctrlPts.std]';
    end
    cpId=prob.ctrlPts(:,1);
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
    %s0.estIO([1:5,7:8],:)=true;
    % Estimate px,py,c,K1-K3,P1-P2.
    s0.estIO([1:3],:)=true;
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
