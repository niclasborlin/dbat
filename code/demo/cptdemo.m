% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

dampings={'none','gna','lm','lmp'};

dampings=dampings(2);

if ~exist('fName','var')
    fName=fullfile(curDir,'data','weighted','w0cm','w0cm-pmexport.txt');
    cpName=fullfile(curDir,'data','weighted','w0cm','ctrlpts.txt');
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
    [cpID,CP,CPS,cpNames]=loadcpt(cpName);
    if ~isequal(sort(cpID(:)),prob.ctrlPts(:,1))
        error('Control point id mismatch.');
    end
    % Replace XYZ positions.
    [~,ia,ib]=intersect(cpID,prob.ctrlPts(:,1));
    prob.ctrlPts(ib,2:4)=CP(:,ia)';
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

fprintf(['Using damping %s. To use another damping, modify line 6 ' ...
         'of %s\n'],dampings{1},mfilename);

% Optionally, change estimate/prior handling.

s0=ss0;
% Reset all parameters to be estimated.
s0.IO(s0.estIO)=nan;
s0.IO(s0.useIOobs)=s0.prior.IO(s0.useIOobs);

s0.EO(s0.estEO)=nan;
s0.EO(s0.useEOobs)=s0.prior.EO(s0.useEOobs);

s0.OP(s0.estOP)=nan;
s0.OP(s0.useOPobs)=s0.prior.OP(s0.useOPobs);

% Use sigma0=.1 pixels as first approximation.
s0.markStd(:)=.1;

s1=resect(s0,'all',cpID,1,0,cpID);
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
        fprintf('Bundle ok after %d iterations with sigma0=%.2f pixels\n', ...
                iters(i),sigma0(i));
    else
        fprintf(['Bundle failed after %d iterations. Last sigma0 estimate=%.2f ' ...
                 'pixels\n'],iters(i),sigma0(i));
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

% EOcmp [Omega, Phi, Kappa, Xc, Yc, Zc]
PMEOvalues=[33.115919
0.774225
-0.807417
0.533043
-0.481070
1.297253
27.519666
1.519317
89.410506
0.548134
-0.588805
1.472052
28.931491
-27.294733
-41.720325
-0.368449
-0.343896
1.254891
28.497681
-23.574097
48.195345
-0.448021
-0.527640
1.613202
-0.541059
-34.010433
-90.720093
-0.513486
0.501924
1.287181
1.084799
-24.895726
0.218419
-0.511651
0.480986
1.568532
-30.232741
-25.470412
-142.481397
-0.358909
1.415510
1.276305
-26.133597
-25.111303
-53.379823
-0.543542
1.474817
1.651559
-34.860522
0.606442
179.922177
0.520274
1.511275
1.275195
-27.440697
-1.857735
-90.754388
0.434815
1.578635
1.467405
-27.916172
27.615928
138.616728
1.338261
1.359611
1.257242
-28.220799
23.209787
-130.936209
1.415573
1.497235
1.610388
-0.695955
34.766836
89.973816
1.535774
0.515393
1.283751
-1.623461
28.277566
180.319791
1.592166
0.540604
1.433149
30.571822
25.215426
36.144163
1.361393
-0.355302
1.305785
26.771485
26.230352
126.504900
1.550245
-0.449245
1.677349
0.877786
3.063271
90.483996
0.558099
0.099209
1.878746
1.219504
3.671129
88.475219
0.653684
0.816553
1.895275
2.399268
2.519687
-0.916526
0.908012
0.465390
1.761481
0.359850
6.962005
-0.642139
0.466417
0.477806
1.794249];

PMEOstd=[
    0.014
0.014
0.008
4.1e-004
3.2e-004
3.0e-004
0.016
0.016
0.008
5.1e-004
4.0e-004
4.4e-004
0.016
0.014
0.009
3.6e-004
3.6e-004
3.9e-004
0.018
0.017
0.010
5.4e-004
5.2e-004
5.8e-004
0.017
0.013
0.008
3.1e-004
4.1e-004
3.0e-004
0.019
0.017
0.008
4.6e-004
5.6e-004
4.5e-004
0.015
0.014
0.009
3.7e-004
3.6e-004
4.0e-004
0.019
0.017
0.010
5.4e-004
5.5e-004
6.0e-004
0.013
0.014
0.008
4.1e-004
3.1e-004
3.0e-004
0.016
0.016
0.008
5.0e-004
4.0e-004
4.3e-004
0.016
0.014
0.009
3.6e-004
3.6e-004
3.9e-004
0.019
0.017
0.010
5.3e-004
5.3e-004
5.7e-004
0.017
0.013
0.008
3.1e-004
4.1e-004
3.1e-004
0.018
0.016
0.008
3.8e-004
4.9e-004
4.3e-004
0.016
0.014
0.009
3.7e-004
3.8e-004
4.0e-004
0.019
0.017
0.010
5.5e-004
5.7e-004
6.1e-004
0.024
0.026
0.007
8.9e-004
8.2e-004
3.3e-004
0.023
0.023
0.007
8.0e-004
7.9e-004
2.8e-004
0.024
0.022
0.007
7.0e-004
7.7e-004
3.0e-004
0.027
0.024
0.006
8.0e-004
8.7e-004
2.0e-004
];

PMEOvalues=reshape(PMEOvalues,6,[]);
PMEOstd=reshape(PMEOstd,6,[]);

D=diag(kron([1,pi/180],ones(1,3)));

PMEOvalues=D*PMEOvalues([4:6,1:3],:);
PMEOstd=D*PMEOstd([4:6,1:3],:);

adv=abs((PMEOvalues-result{i}.EO(1:6,:)));
adv=min(adv,abs(adv-2*pi));
maxEOvalDiff=max(adv(:))

EOstd=sqrt(full(reshape(diag(bundle_cov(result{1},E{1},'CEO')),6,[])));
EOstdDiff=abs(PMEOstd-EOstd);
maxEOstdDiff=max(EOstdDiff(:))

