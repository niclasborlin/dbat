disp(['Bundle adjustment test, real data from loaded photomodeler' ...
	  ' export file'])

try
    % Provoke early error if ptb not in path.
    PTGaussian(1);
catch
    addpath(fullfile(niclas,'ptb','code','ptb'),'-end');
end

%dataSets=repmat({1:60,1:4:60,1:20,3:17,1:2:60},1,4);
dataSets=repmat({1:60,1:4:60,1:20,1:2:20,1:5},1,4);

trim5degpts=repmat([false,false,false,false,false,true,true,true,true,true],1,2);
freeNet=[false(1,10),true(1,10)];

seed=24;

rand('seed',seed);

%dataSetNo=12
%veto=true
dataSetNo
changeBadPts
removeBadPts
veto
vetoAll
largePert
mirror
doSave
angThres=5*pi/180;
vetoTiming
closeLimit1

doPlot
doPlotX0veto
reorder=false
calcProblemData

freeNetAdjustment=freeNet(dataSetNo)

if veto
    fileName=['romabundle_set',int2str(dataSetNo),'_veto'];
    if vetoAll
        fileName=[fileName,'_all'];
    end
    if mirror
        fileName=[fileName,'_mirror'];
    end
else
    fileName=['romabundle_set',int2str(dataSetNo)];
end

if largePert
    fileName=[fileName,'_large'];
end

if smallPert
    fileName=[fileName,'_small'];
end

if removeBadPts
    fileName=[fileName,'_removebad'];
end

if changeBadPts
    fileName=[fileName,'_changebad'];
end

host=getenv('HOST');

D=22; % Object diameter.

if smallPert
    % Position perturbations, percent.
    cPert=[0:5]/100;
    % Angular perturbations, degrees.
    angPert=[0:0.5:3]/180*pi;
else
    % Position perturbations, percent.
    cPert=[0:5]/100;
    % Angular perturbations, degrees.
    angPert=[0:5]/180*pi;
end
% Simulations per level.
nRuns=250;
if doBench
    cIx=2;
    aIx=1;
    cPert=cPert(cIx);
    angPert=angPert(aIx);
    fileName=[fileName,'_bench_',host];
end

fileName=[fileName,'.mat'];

if largePert
    cPert=cPert*2;
    %angPert=angPert*2;
end

dataset=dataSets{dataSetNo};

if exist(fileName,'file') && ~calcProblemData
    Z=load(fileName);
    savedResults=Z.results;
    if ~isfield(savedResults,'host')
        savedResults.host='unknown';
    end
    if ~isfield(savedResults,'seed')
        savedResults.seed=24;
    end
    % Calculate # of completed iterations.
    sz=size(savedResults.time);
    if length(sz)<4
        sz(end+1:4)=1;
    end
    t=permute(savedResults.time,[1,3,4,2]);
    lastIter=min([find(any(isnan(reshape(t,prod(sz([1,3,4])),[])))),inf])-1;
    if isinf(lastIter)
        lastIter=size(savedResults.time,2);
    end
    nIter=sz(2);
    fprintf('Loaded data from %s: ',fileName);
    fprintf('D = %g, nRuns=%d, %d done (%d%%), host=%s.\n',savedResults.D,...
            nIter,lastIter,floor(lastIter/nIter*100),savedResults.host);
    if ~isequal(D,savedResults.D)
        warning('Different D')
    end
    if ~isequal(cPert,savedResults.cPert)
        warning('Different cPert')
    end
    if ~isequal(angPert,savedResults.angPert)
        warning('Different angPert')
    end
    if ~isequal(veto,savedResults.veto)
        warning('Different veto')
    end
    if ~isequal(dataset,savedResults.dataset)
        warning('Different dataset')
    end
    if ~isequal(seed,savedResults.seed)
        warning('Different seed')
    end
    if ~isequal(host,savedResults.host)
        warning('Different host')
    end
    if doBench
        yn='y';
    else
        yn=input('Do you want to use stored data (y/n)? ','s');
    end
    switch yn
    case 'y'
        % Everything is already done.
    case 'n'
        savedResults=[];
        lastIter=0;
    otherwise
        error('Bad reply')
    end
else
    savedResults=[];
    lastIter=0;
end

if ~exist('prob','var')
    prob=loadpm('../data/roma.txt',true);
    drawnow
end

if all(~isnan(prob.imSz))
    imSz=prob.imSz;
else
    error('Unknown image size');
end

% Construct internal parameters.
f=prob.cameras(1).inner(1);
wh=prob.cameras(1).inner(4:5);
pp=prob.cameras(1).inner(2:3);
pp(2)=-pp(2);
K=-prob.cameras(1).inner(6:7);
P=-prob.cameras(1).inner(9:10);
a=zeros(2,1);
u=repmat(mean(imSz(:)./wh(:)),2,1);

nK=length(K);
nP=length(P);

IO=[pp(:);f;K(:);P(:);a(:);u];

camCalMat=diag([-f,-f,1]);
camCalMat(1:2,3)=pp';

useCams=dataset;

fixCam=min(useCams);

% Data storage:
%
% IO internal parameters.
% EO external parameters, columnwise.
%
% OPid contains object point ids
% OP   contains object point XYZ, sorted by id.
%
% isCtrl shows if OP is control point
%
% vis(i,j)=1 if object point i is visible in image j
%
% markPos mark positions sorted by image, id.


% Camera stations.
cams=cat(1,prob.cameras.outer);
C=cams(:,1:3)';
ang=cams(:,[6,5,4])'/180*pi;

% Construct exterior orientation parameters.
EO=[C;ang];
% Indicate all rotations are Euler 'xyz'.
EO(7,1)=0;

% Object points.
OP=prob.objPts(:,2:4)';
OPid=prob.objPts(:,1);

% Transform EO, OP to have EO(fixCam) defining the origin.

% Shift fixCam to origin.
sceneC=EO(1:3,fixCam);
EO(1:3,:)=EO(1:3,:)-repmat(sceneC,1,size(EO,2));
OP=OP-repmat(sceneC,1,size(OP,2));

% Rotate 
RR=pm_eulerrotmat(EO(4:6,fixCam));
%RR=eye(3);
% Rotate OP
OP=RR*OP;
% Rotate camera positions
EO(1:3,:)=RR*EO(1:3,:);
for i=1:size(EO,2)
    R=pm_eulerrotmat(EO(4:6,i));
    R2=R*RR';
    ang=derotmat3d(R2);
    EO(4:6,i)=ang(:);
end

if all(abs(EO(4:6,fixCam)<eps*10))
    EO(4:6,fixCam)=0;
end
    
% Sort according to id.
[dummy,ix]=sort(OPid);
OPid=OPid(ix);
OP=OP(:,ix);
if (isempty(prob.ctrlPts))
	isCtrl=false(size(OPid));
else
	isCtrl=ismember(OPid,prob.ctrlPts(:,1));
end

% Sort mark points after image, then id.
markPts=msort(prob.markPts);

% Remove mark points in unused images.
imNo=markPts(:,1)'+1;
unused=~ismember(imNo,useCams);
markPts(unused,:)=[];

% Remove unreconstructed points.
markId=markPts(:,2);
reconstructed=ismember(markId,OPid);
markPts=markPts(reconstructed,:);

% Remove points visible in too few images.
imNo=markPts(:,1)'+1;
markId=markPts(:,2);
vis=sparse(markId,imNo,true,max(OPid),max(useCams)); % Mapping from (id,image)
vis=vis(OPid,:); % Mapping from (OPid(i),image)
% Keep OP visible in more than 1 image and any visible ctrl pts.
ptOK=sum(vis,2)>1 | sum(vis,2)>0 & isCtrl;

% Optionally, remove all points <5 deg intersection angle.
if trim5degpts(dataSetNo)
    maxAng=-inf(size(isCtrl));
    visT=vis';
    % For each row of vis...
    for i=1:size(vis,1)
        j=find(visT(:,i));
        % We need at least two observations...
        if length(j)>1 && ~isCtrl(i)
            % Calculate intersection angle.
            dpt=repmat(OP(:,i),1,length(j))-EO(1:3,visT(:,i));
            for k1=1:length(j)-1
                for k2=k1+1:length(j)
                    ang=subspace(dpt(:,k1),dpt(:,k2));
                    maxAng(i)=max(maxAng(i),ang);
                end
            end
        end
    end
    angOK=maxAng>angThres | isCtrl;
    min(maxAng(isfinite(maxAng)))*180/pi
    min(maxAng(maxAng>angThres))*180/pi
else
    angOK=true(size(isCtrl));
end

ptOK=ptOK & angOK;

OPid=OPid(ptOK);
OP=OP(:,ptOK);
isCtrl=isCtrl(ptOK);
vis=vis(ptOK,:);

% Only keep marks that belong to ok OP.
markId=markPts(:,2);
markPtOK=ismember(markId,OPid);
markPts=markPts(markPtOK,:);

imNo=markPts(:,1)'+1;
markId=markPts(:,2);
markPos=markPts(:,3:4)';

% Create mapping from OPid,image to column in markPos.
colPos=reshape(cumsum(vis(:)),size(vis)).*vis;

% Number of object points.
n=size(OP,2);

% Number of cameras.
m=size(EO,2);

disp('Using nominal -y')
markPos(2,:)=-markPos(2,:);

correctedPt=roll(pm_multilenscorr1(markPos,IO,nK,nP,1,1),2);

OPsave=OP;
EOsave=EO;
IOsave=IO;

% Don't estimate any IO.
cIO=false(size(IO));

% Estimate all OP except CP.
cOP=true(size(OP));
cOP(:,isCtrl)=0;

% Estimate only used cameras...
cEO=false(size(EO));
cEO(1:6,useCams)=true;

% ...except the fixed camera....
cEO(:,fixCam)=false;
% ...and the largest coordinate shift among the other cameras...
d=EO(1:3,:)-repmat(EO(1:3,fixCam),1,size(EO,2));
[dummy,i,j]=max2(abs(d.*cEO(1:3,:)));
cEO(i,j)=false;

OP=pm_multiforwintersect(IO,EO,1,colPos,correctedPt,find(any(cOP,1)));

intersectAngle=nan(1,size(cOP,2));
% Pre-transpose for speed.
visT=vis';
for i=find(any(cOP,1))
    % For each row of vis...
    j=find(visT(:,i));
    % We need at least two observations...
    if length(j)>1
        % Calculate intersection angle.
        maxAng=-inf;
        dpt=repmat(OP(:,i),1,length(j))-EO(1:3,visT(:,i));
        for k1=1:length(j)-1
            for k2=k1+1:length(j)
                ang=subspace(dpt(:,k1),dpt(:,k2));
                maxAng=max(maxAng,ang);
            end
        end
        intersectAngle(i)=maxAng;
    end
end

d=pm_multidepth(IO,EO,OP,vis,1);

PI=zeros(4,0);
constr=[];
orthoP=zeros(0,2);
cPI=[];
cSI=[];
oIO=[];
oEO=[];
oOP=[];
vIO=[];
vEO=[];
vOP=[];
SI=[];
shiftVal=[];
shiftIx=[];
ctrlPtFix=[];

base=1;
if ~reorder
    % Original ordering
    [ixIO,base]=pindex(nnz(cIO),base);
    [ixEO,base]=pindex(nnz(cEO),base);
    [ixOP,base]=pindex(nnz(cOP),base);
    [ixPI,base]=pindex(nnz(cPI),base);
    [ixSI,base]=pindex(nnz(cSI),base);
    fun='pm_eulerbundle1p';
else
    % Reordered.
    [ixSI,base]=pindex(nnz(cSI),base);
    [ixPI,base]=pindex(nnz(cPI),base);
    [ixOP,base]=pindex(nnz(cOP),base);
    [ixEO,base]=pindex(nnz(cEO),base);
    [ixIO,base]=pindex(nnz(cIO),base);
    fun='pm_eulerbundle1r';
end

% Run one bundle to get true values.
x0=zeros(nnz(cSI)+nnz(cPI)+nnz(cOP)+nnz(cEO)+nnz(cIO),1);
x0(ixIO)=IO(cIO);
x0(ixEO)=EO(cEO);
x0(ixOP)=OP(cOP);
x0(ixPI)=PI(cPI);
x0(ixSI)=SI(cSI);

params={markPos,1,IO,nK,nP,EO,1,OP,vis,PI,constr,orthoP,cIO,cEO,cOP,cPI,cSI,oIO,oEO,oOP,vIO,vEO,vOP,shiftVal,shiftIx,ctrlPtFix};

if (1)
	disp('Damped solution');
	mu=1e-1;
else
	disp('Undamped solution')
	mu=-inf;
end

mu=0.1;
maxIter=20;
W=speye(prod(size(markPos))+nnz(vIO)+nnz(vEO)+nnz(vOP));
tolf=1e-3;
alphaMin=1e-3;
tolc=1e-6;
nu0=0.1;
stopWatch=cputime;
[x,nIter,code,lambda,X,alphas,cc,ll,nus,F,J]=sqp_paper(fun,'',x0,nu0,mu,maxIter,tolf,tolc,W,params);
time=cputime-stopWatch
code,nIter,alphas

asfds
% [f,J]=pm_eulerbundle1p_f(x,params{:});
%K=pm_eulerbundle1_n(x,params{:});
%k=size(K,1);
sigma0=sqrt((F'*F)/(size(J,1)-size(J,2)));
sigma0px=sigma0*mean(u)
if calcProblemData
    Qxxpart=sigma0^2*((J'*J)\speye(size(J,2),nnz(cEO)));
    EOC=PTGaussian(zeros(size(EO)));
    EOC(cEO)=PTGaussian(x(1:nnz(cEO)),Qxxpart(1:nnz(cEO),1:nnz(cEO)));
end

T0=blkdiag(1,[0,-1;1,0],1);
camSize=D/22*10;

if doPlot || calcProblemData
    fig=tagfigure('romabundle');
    [EOT,OPT]=pm_multixform(EO,OP,blkdiag(pm_eulerrotmat([-15,0,0]*pi/180),1));
    pm_plotmulti(IO,EOT,OPT,isCtrl,any(vis,1),cIO,cEO,cOP,[],ixIO,ixEO,ixOP,{}, ...
                 camSize,fig,T0);
    drawnow
    set(findobj(fig,'type','line'),'markersize',0.5);
    axis off
    title('');
    view(-30,14);

    if max(dataset)>50
        dataset=dataset([1:end,1]);
    end
    trueAnglesLeft=nan(3,max(dataset));
    trueAnglesRight=nan(3,max(dataset));
    for i=2:length(dataset)
        cref=dataset(i-1);
        ci=dataset(i);
        [EOA,OPA]=pm_multialign(EO,OP,cref);
        trueAnglesRight(:,ci)=EOA(4:6,ci);
        [EOA,OPA]=pm_multialign(EO,OP,ci);
        trueAnglesLeft(:,cref)=EOA(4:6,cref);
    end
    
    anglesLeft=nan(3,max(dataset));
    anglesRight=nan(3,max(dataset));
    % Determine camera angles based on Nister's 5-point algorithm.
    for ii=2:length(dataset)
        ii
        cam1=dataset(ii-1);
        cam2=dataset(ii);
        % Points visible in these cameras.
        i=all(vis(:,[cam1,cam2]),2);
        % Extract mm coordinates.
        pts1=correctedPt(:,colPos(i,cam1));
        pts2=correctedPt(:,colPos(i,cam2));
        % Normalize
        pt1n=camCalMat\homogeneous(pts1);
        pt2n=camCalMat\homogeneous(pts2);
        Evec=essmat5(pt2n,pt1n);
        bestSp=-inf;
        bestErr=inf;
        PP1=cell(1,size(Evec,2));
        PP2=cell(1,size(Evec,2));
        XXX=cell(1,size(Evec,2));
        IIF=cell(1,size(Evec,2));
        SP=nan(4,size(Evec,2));
        for i=1:size(Evec,2)
            [P1,P2,Xnew,inFront,sp]=camsfrome(reshape(Evec(:,i),3,3),...
                                              pt1n,pt2n,-1);
            PP1{i}=P1;
            PP2{i}=P2;
            XXX{i}=Xnew;
            IIF{i}=inFront;
            SP(:,i)=sp;
            if 0
                C2=euclidean(null(P2));
                R2=P2(:,1:3);
                ang=derotmat3d(R2)';
                EEO=[zeros(7,1),[C2;derotmat3d(R2)';0]];
                ff=figure(tagfigure('nister'));
                pm_plotmulti([],EEO,Xnew(1:3,:),false(1,size(Xnew,2)),true(1,2),...
                             [],[],[],[],[],[],[],{},camSize,ff,T0);
            end
        end
        D=SP(1,:)-SP(2,:);
        best=find(D==max(D));
        if SP(2,best)>0
            best=[];
        elseif length(best)>1
            bestErr=zeros(size(best));
            bestErrMax=zeros(size(best));
            for i=1:length(best)
                ci=best(i);
                P1=PP1{ci};
                P2=PP2{ci};
                Xnew=homogeneous(XXX{ci});
                bestErr(i)=sum(sum((euclidean(P1*Xnew)-euclidean(pt1n)).^2))+...
                    sum(sum((euclidean(P2*Xnew)-euclidean(pt2n)).^2));
                bestErrMax(i)=max(...
                    max(sum((euclidean(P1*Xnew)-euclidean(pt1n)).^2)),...
                    max(sum((euclidean(P2*Xnew)-euclidean(pt2n)).^2)));
                
            end
            [dummy,i]=min(bestErr);
            best=best(i(1));
        end
        if ~isempty(best)
            P2=PP2{best};
            C2=euclidean(null(P2));
            R2=P2(:,1:3);
            anglesRight(:,cam2)=derotmat3d(R2)';
            EEO=[zeros(7,1),[C2;anglesRight(:,cam2);0]];
            EEOA=pm_multialign(EEO,[],2);
            anglesLeft(:,cam1)=EEOA(4:6,1);
        end
        (anglesRight-trueAnglesRight)*180/pi
    end
    figure(tagfigure('plots'))
    subplot(2,1,1)
    hist(intersectAngle*180/pi,1:90)
    title('Intersection angles')
    subplot(2,1,2)
    hist(sum(vis,2),2:max(sum(vis,2)))
    title('Ray count')
    disp(sprintf('Num cameras: %d',nnz(dataset)));
    disp(sprintf('Num OP:      %d',size(OP,2)));
    disp(sprintf('Avg ray cnt: %.1f',full(mean(sum(vis,2)))));
    disp(sprintf('2 rays:      %.0f%%',nnz(sum(vis,2)==2)/size(OP,2)*100));
    disp(sprintf('avg alpha:   %.0f',mean(intersectAngle)*180/pi));
    disp(sprintf('min alpha:   %.2g',min(intersectAngle)*180/pi));
    disp(sprintf('# alpha<5:   %.0f%%',nnz(intersectAngle<angThres)/size(OP,2)*100));
    disp(sprintf('Sigma0:      %.2f',sigma0px));
end

% Store true values.
IOtrue=IO;
EOtrue=EO;
OPtrue=OP;
PItrue=PI;
SItrue=SI;
IOtrue(cIO)=x(ixIO);
EOtrue(cEO)=x(ixEO);
OPtrue(cOP)=x(ixOP);
PItrue(cPI)=x(ixPI);
SItrue(cSI)=x(ixSI);
xTrue=x;

if veto
    if vetoAll
        vetoFun='pm_eulerbundle1p_v3all';
    else
        vetoFun='pm_eulerbundle1p_v3';
    end
else
    vetoFun='';
end
vetoRange=[0,inf];
%vetoX0Range=[1,2*D];
vetoX0Range=vetoRange;
if closeLimit1
    vetoX0Range=[1,inf];
    fileName=[fileName,'_close1'];
end

if calcProblemData
    return;
end
    
methodsToTest=1:4;

nMethods=length(methodsToTest);
ncPert=length(cPert);
nAngPert=length(angPert);

total=nRuns*ncPert*nAngPert;
done=0;
startTime=now;

results=struct('iters',zeros(nMethods,nRuns,ncPert,nAngPert,'int8'),...
               'code',zeros(nMethods,nRuns,ncPert,nAngPert,'int8'),...
               'time',nan(nMethods,nRuns,ncPert,nAngPert,'single'),...
               'diffNormAvg',nan(nMethods,nRuns,ncPert,nAngPert,'single'),...
               'diffNormMax',nan(nMethods,nRuns,ncPert,nAngPert,'single'),...
               'host',host,'seed',seed,...
               'badPts',zeros(nMethods,nRuns,ncPert,nAngPert,'uint32'),...
               'D',D,'nRuns',nRuns,'total',total,'done',done,'veto',veto,...
               'vetox0',zeros(nMethods,nRuns,ncPert,nAngPert,'uint32'),...
               'vetoActivated',false(nMethods,nRuns,ncPert,nAngPert),...
               'cPert',cPert,'angPert',angPert,'dataset',dataset);

if ~isempty(savedResults)
    % # of runs in previously stored data.
    preRuns=size(savedResults.iters,2);
    mm=size(savedResults.iters,3);
    nn=size(savedResults.iters,4);
    results.iters(:,1:preRuns,1:mm,1:nn)=savedResults.iters;
    results.code(:,1:preRuns,1:mm,1:nn)=savedResults.code;
    results.time(:,1:preRuns,1:mm,1:nn)=savedResults.time;
    results.vetox0(:,1:size(savedResults.vetox0,2),1:mm,1:nn)=savedResults.vetox0;
    results.diffNormAvg(:,1:preRuns,1:mm,1:nn)=savedResults.diffNormAvg;
    results.diffNormMax(:,1:preRuns,1:mm,1:nn)=savedResults.diffNormMax;
    results.D=savedResults.D;
    results.nRuns=nRuns;
    results.total=total;
    if ~strcmp(results.host,host)
        results.host='mixed';
    end
    preDone=savedResults.done;
else
    preDone=0;
end

done=preDone;

% Find out which iterations we haven't done.
[i1,i2,i3,i4]=ind2sub(size(results.time),find(isnan(results.time)));

firstIter=min(i2);

% Unscaled perturbations.
pert=rand(6,size(EO,2),nRuns);

lastMethod=0;

if freeNetAdjustment
    % No fixed EO parameters in free net adjustment.
    cEO(1:6,useCams)=true;
    base=1;
    [ixIO,base]=pindex(nnz(cIO),base);
    [ixEO,base]=pindex(nnz(cEO),base);
    [ixOP,base]=pindex(nnz(cOP),base);
    [ixPI,base]=pindex(nnz(cPI),base);
    [ixSI,base]=pindex(nnz(cSI),base);
    params={markPos,1,IO,nK,nP,EO,1,OP,vis,PI,constr,orthoP,cIO,cEO,cOP,cPI,cSI,oIO,oEO,oOP,vIO,vEO,vOP,shiftVal,shiftIx,ctrlPtFix};
end

cOPglobal=cOP;
visGlobal=vis;
markPosGlobal=markPos;
ixTrue=1:length(x);
ixOPglobal=ixOP;

for ii=firstIter:nRuns
    for cPertLevel=ncPert:-1:1
        for angPertLevel=nAngPert:-1:1
            if ~isnan(results.time(1,ii,cPertLevel,angPertLevel))
                fprintf('Already done iter %d, cLevel %d, angLevel %d.\n',...
                        ii,cPertLevel,angPertLevel);
                continue
            end
            IOtry=IOtrue;
            EOtry=EOtrue;
            OPtry=OPtrue;
            PItry=PItrue;
            SItry=SItrue;
            cOP=cOPglobal;
            ixOP=ixOPglobal;
            vis=visGlobal;
            markPos=markPosGlobal;

            % Create perturbation.
            EOp=[D*cPert(cPertLevel)*(2*pert(1:3,:,ii)-1);
                angPert(angPertLevel)*(2*pert(4:6,:,ii)-1)];
            
            %IOtry(cIO)=
            EOp(end+1,1)=0;
            EOtry(cEO)=EOtry(cEO)+EOp(cEO);
    
            OPtry=pm_multiforwintersect(IO,EOtry,1,colPos,correctedPt,...
                                        find(any(cOP,1)));
            OPNoVetoX0=OPtry;

            if changeBadPts || removeBadPts
                d=pm_multidepth(IO,EOtry,OPtry,vis,1);
                bad=d<0;
            end

            ixCmpTrue=ixTrue;
            if removeBadPts
                if nnz(bad)>0
                    % Ignore bad OP
                    badOP=any(bad,2);
                    cOP(:,badOP)=false;
                    ixOP3=reshape(ixOPglobal,3,[]);
                    ixCmpTrue=setdiff(ixCmpTrue,reshape(ixOP3(:,badOP),1,[]));
                    base=1;
                    [ixIO,base]=pindex(nnz(cIO),base);
                    [ixEO,base]=pindex(nnz(cEO),base);
                    [ixOP,base]=pindex(nnz(cOP),base);
                    [ixPI,base]=pindex(nnz(cPI),base);
                    [ixSI,base]=pindex(nnz(cSI),base);
                    % New visibility matrix.
                    vis=vis & ~bad;
                    % Remove ignored observations.
                    keep=colPos.*~bad;
                    markPos=markPos(:,keep(keep~=0));
                end
            end
            W=speye(prod(size(markPos))+nnz(vIO)+nnz(vEO)+nnz(vOP));
            params={markPos,1,IO,nK,nP,EO,1,OP,vis,PI,constr,orthoP,cIO,cEO,cOP,cPI,cSI,oIO,oEO,oOP,vIO,vEO,vOP,shiftVal,shiftIx,ctrlPtFix};
            
            if changeBadPts
                if doPlotX0veto
                    vFig=tagfigure('vetox0');
                    wFig=tagfigure('vetox0pt');
                end                    
                if vetoAll
                    visX0=repmat(any(vis,1),size(vis,1),1);
                else
                    visX0=vis;
                end
                for repeat=1:5
                    if mirror
                        d=pm_multidepth(IO,EOtry,OPtry,vis,1);
                        bad=d<0;
                    else
                        d=pm_multidepth(IO,EOtry,OPtry,visX0,1);
                        bad=d<vetoX0Range(1) | d>vetoX0Range(2);
                    end
                    
                    if doPlotX0veto
                        pm_plotmulti(IO,EOtry,OPtry,isCtrl,any(vis,1),...
                                     cIO,cEO,cOP,[],ixIO,ixEO,ixOP,{},...
                                     camSize,vFig,T0);
                        str=sprintf('%d bad points in cameras',...
                                    nnz(any(bad,2)));
                        str=[str,sprintf(' %d',find(any(bad,1)))];
                        title(gca(vFig),str);
                        pause
                    end
                    if ~any(bad(:))
                        break;
                    else
                        if mirror
                            % For each point behind...
                            for i=find(any(bad,2))'
                                % Find out which cameras the point is behind...
                                cams=find(bad(i,:));
                                % Center of gravity of those cameras...
                                cg=mean(EOtry(1:3,cams),2);
                                % Reflect point.
                                oldPt=OPtry(:,i);
                                newPt=cg+(cg-OPtry(:,i));
                                OPtry(:,i)=newPt;
                                if false && doPlotX0veto
                                    OPP=[oldPt,newPt,EOtry(1:3,vis(i,:))];
                                    l={};
                                    for j=3:size(OPP,2)
                                        l{end+1}=[1,j];
                                        l{end+1}=[2,j];
                                    end
                                    pm_plotmulti(IO,EOtry,OPP+eps,...
                                                 false(1,size(OPP,2)),...
                                                 vis(i,:),cIO,cEO,...
                                                 cOP,[],ixIO,ixEO,ixOP,l,...
                                                 camSize,wFig,T0);
                                    str=sprintf('id=%d, cams=',i);
                                    str=[str,sprintf('%d,',find(vis(i,:)))];
                                    str=[str,' d_b='];
                                    str=[str,sprintf('%.2f,',...
                                                     full(d(i,vis(i,:))))];
                                    dNew=pm_multidepth(IO,EOtry,newPt,...
                                                       vis(i,:),1);
                                    str=[str,' d_a='];
                                    str=[str,sprintf('%.2f,',...
                                                     full(dNew(vis(i,:))))];
                                    title(gca(wFig),str);
                                    pause
                                end
                            end
                        else
                            % For each bad point...
                            for i=find(any(bad,2))'
                                % Find out the bad cameras w.r.t. this point
                                cams=find(vis(i,:));
                                if nnz(cams)~=nnz(find(bad(i,:)))
                                    %disp('Veto x0 diff');
                                end
                                % Find all good points visible in the
                                % same cameras.
                                ptGood=all(vis(:,cams) & ~bad(:,cams),2);
                                % Set bad point to average of good points.
                                if nnz(ptGood)==0
                                    warning('BAD');
                                end
                                OPtry(:,i)=mean(OPtry(:,ptGood),2);
                            end
                        end
                    end
                end

                if mirror
                    for repeat=1:5
                        d=pm_multidepth(IO,EOtry,OPtry,visX0,1);
                        bad=d<vetoX0Range(1);
                    
                        if doPlotX0veto
                            pm_plotmulti(IO,EOtry,OPtry,isCtrl,any(vis,1),...
                                         cIO,cEO,cOP,[],ixIO,ixEO,ixOP,{},...
                                         camSize,vFig,T0);
                            str=sprintf('%d bad points in cameras',...
                                        nnz(any(bad,2)));
                            str=[str,sprintf(' %d',find(any(bad,1)))];
                            title(gca(vFig),str);
                            pause
                        end
                        if ~any(bad(:))
                            break;
                        end
                        % For each point too close to a camera...
                        for i=find(any(bad,2))'
                            % Find out which cameras the point is behind...
                            cams=find(bad(i,:));
                            for j=cams
                                % Calculate direction of principal axis.
                                R=rotmat(EO(4:6,j));
                                pa=R(3,:)';
                                pas=pa*(vetoX0Range(1)-d(i,j));
                                % Push point forward.
                                oldPt=OPtry(:,i);
                                newPt=OPtry(:,i)-pas;
                                OPtry(:,i)=newPt;
                                if false && doPlotX0veto
                                    OPP=[oldPt,newPt,...
                                         EOtry(1:3,j),EOtry(1:3,j)-pa,...
                                         EOtry(1:3,vis(i,:) | bad(i,:))];
                                    l={[3,4]};
                                    for j=5:size(OPP,2)
                                        l{end+1}=[1,j];
                                        l{end+1}=[2,j];
                                    end
                                    pm_plotmulti(IO,EOtry,OPP+eps,...
                                                 false(1,size(OPP,2)),...
                                                 vis(i,:) | bad(i,:),cIO,cEO,...
                                                 cOP,[],ixIO,ixEO,ixOP,l,...
                                                 camSize,wFig,T0);
                                    str=sprintf('id=%d, cams=',i);
                                    str=[str,sprintf('%d,',...
                                                     find(vis(i,:) | bad(i,:)))];
                                    str=[str,' d_b='];
                                    str=[str,sprintf('%.2f,',...
                                                     full(d(i,vis(i,:)|bad(i,:))))];
                                    dNew=pm_multidepth(IO,EOtry,newPt,...
                                                       vis(i,:)|bad(i,:),1);
                                    str=[str,' d_a='];
                                    str=[str,sprintf('%.2f,',...
                                                     full(dNew(vis(i,:)|bad(i,:))))];
                                    title(gca(wFig),str);
                                    pause
                                end
                            end
                        end
                    end 
                end
            end
            
            % How many points were affected by the veto X0 algorithm?
            nVetoX0=nnz(any(OPtry~=OPNoVetoX0));
            
            % If free net adjustment, shift center of object to origin.
            if freeNetAdjustment
                c=mean(OPtry,2);
                OPtry=OPtry-repmat(c,1,size(OPtry,2));
                EOtry(1:3,:)=EOtry(1:3,:)-repmat(c,1,size(EOtry,2));
            end
            
            x0=zeros(nnz(cSI)+nnz(cPI)+nnz(cOP)+nnz(cEO)+nnz(cIO),1);
            x0(ixIO)=IOtry(cIO);
            x0(ixEO)=EOtry(cEO);
            x0(ixOP)=OPtry(cOP);
            x0(ixPI)=PItry(cPI);
            x0(ixSI)=SItry(cSI);
            
            for mi=[2,1,3:length(methodsToTest)]
                mm=methodsToTest(mi);
                if veto
                    feval(vetoFun,'setup',vetoRange,{},[],[],OPid,dataset);
                end
                switch mm
                case 1
                    if lastMethod==2 && all(alphas==1) && ~vetoTiming
                        % Do nothing, results of damped will be same as
                        % undamped, unless vetoTiming is requested.
                    else
                        mu=-inf;
                        if freeNet(dataSetNo)
                            stopWatch=cputime;
                            [x,code,nIter,X,alphas,F]=...
                                gaussn_paper_fn('pm_eulerbundle1p_f',...
                                                'pm_eulerbundle1p_fnc',...
                                                '',x0,...
                                                tolf,maxIter,alphaMin,mu,W,...
                                                params);
                            time=cputime-stopWatch;
                        else
                            stopWatch=cputime;
                            [x,code,nIter,X,alphas,F]=...
                                gaussn_paper('pm_eulerbundle1p_f','',...
                                             x0,tolf,...
                                             maxIter,alphaMin,mu,W,params);
                            time=cputime-stopWatch;
                        end
                    end
                case 2
                    mu=0.1;
                    if freeNet(dataSetNo)
                        stopWatch=cputime;
                        [x,code,nIter,X,alphas,F]=...
                            gaussn_paper_fn('pm_eulerbundle1p_f',...
                                            'pm_eulerbundle1p_fnc','',x0,...
                                            tolf,maxIter,alphaMin,mu,W,...
                                            params);
                        time=cputime-stopWatch;
                    else
                        stopWatch=cputime;
                        [x,code,nIter,X,alphas,F]=...
                            gaussn_paper('pm_eulerbundle1p_f',vetoFun,x0,tolf,...
                                         maxIter,alphaMin,mu,W,params);
                        time=cputime-stopWatch;
                    end
                case 3
                    [f0,J0]=pm_eulerbundle1p_f(x0,params{:});
                    lambda0=1e-10*trace(J0'*J0)/size(J0,2);
                    %lambda0=1e-10;
                    lambdaMin=lambda0;
                    lambdaMax=1e20*lambdaMin;
                    if freeNet(dataSetNo)
                        stopWatch=cputime;
                        [x,code,nIter,X,lambdas,rhos,F]=...
                            levenberg_algebraic_paper_fn(...
                                'pm_eulerbundle1p_f','pm_eulerbundle1p_fnc',vetoFun,x0,tolf,maxIter,...
                                lambda0,lambdaMin,lambdaMax,params);
                        time=cputime-stopWatch;
                    else
                        stopWatch=cputime;
                        [x,code,nIter,X,lambdas,rhos,F]=...
                            levenberg_algebraic_paper(...
                                'pm_eulerbundle1p_f',vetoFun,x0,tolf,maxIter,...
                                lambda0,lambdaMin,lambdaMax,params);
                        time=cputime-stopWatch;
                    end
                case 4
                    delta0=norm(x0);
                    lmmu=0.25;
                    lmeta=0.75;
                    if freeNet(dataSetNo)
                        stopWatch=cputime;
                        [x,code,nIter,X,deltas,rhos,F]=...
                            levenberg_paper_fn('pm_eulerbundle1p_f',...
                                               'pm_eulerbundle1p_fnc',...
                                               vetoFun,x0,...
                                            tolf,maxIter,delta0,lmmu,...
                                            lmeta,params);
                        time=cputime-stopWatch;
                    else
                        stopWatch=cputime;
                        [x,code,nIter,X,deltas,rhos,F]=...
                            levenberg_paper('pm_eulerbundle1p_f',vetoFun,x0,...
                                            tolf,maxIter,delta0,lmmu,...
                                            lmeta,params);
                        time=cputime-stopWatch;
                    end
                end
                lastMethod=mm;
                if code==0 && max(abs(x))>1e+6
                    code=-4;
                end
                errNormAvg=nan;
                errNormMax=nan;
                if freeNetAdjustment
                    % Converged values.
                    EOS=EOtry;
                    OPS=OPtry;
                    EOS(cEO)=x(ixEO);
                    OPS(cOP)=x(ixOP);
                    % Align with camera 1.
                    [EOA,OPA]=pm_multialign(EOS,OPS,1);
                    % Determine scale.
                    s=sqrt(mean(OPtrue(:).^2)/mean(OPA(:).^2));
                    OPA=OPA*s;
                    if code==0
                        errNormAvg=norm(OPtrue(:)-OPA(:))/prod(size(OPA));
                        errNormMax=max(abs(OPtrue(:)-OPA(:)));
                    end
                else
                    if code==0
                        errNormAvg=norm(x-xTrue(ixCmpTrue))/length(x);
                        errNormMax=max(abs(x-xTrue(ixCmpTrue)));
                    end
                end
                if veto
                    vetoActivated=feval(vetoFun,'wasvetoactivated')
                else
                    vetoActivated=false;
                end
                
                results.code(mi,ii,cPertLevel,angPertLevel)=code;
                results.iters(mi,ii,cPertLevel,angPertLevel)=nIter;
                results.time(mi,ii,cPertLevel,angPertLevel)=time;
                results.diffNormAvg(mi,ii,cPertLevel,angPertLevel)=errNormAvg;
                results.diffNormMax(mi,ii,cPertLevel,angPertLevel)=errNormMax;
                results.vetox0(mi,ii,cPertLevel,angPertLevel)=nVetoX0;
                results.vetoActivated(mi,ii,cPertLevel,angPertLevel)=vetoActivated;
                results.badPts(mi,ii,cPertLevel,angPertLevel)=nnz(cOP~=cOPglobal)/3;
                if doPlot
                    if code==0
                        str=sprintf('Method %d, code %d, errAvg=%.2g, errMax=%.2g: ',mm,code,errNormAvg,errNormMax);
                    else
                        str=sprintf('Method %d, code %d: ',mm,code);
                    end
                    pm_plotmulti(IO,EO,OP,isCtrl,any(vis,1),cIO,cEO,cOP,X,...
                                 ixIO,ixEO,ixOP,{},camSize,fig,T0,str);
                    drawnow
                end
            end
            done=done+1;
            t=now;
            etf=(t-startTime)/(done-preDone)*(total-preDone)+startTime;
            fprintf('Data set %d, veto %d, vetoAll %d, mirror %d, largePert %d, removeBadPts=%d, changeBadPts=%d, closeLimit1=%d, host %s.\n',...
                    dataSetNo,veto,vetoAll,mirror,largePert,removeBadPts,changeBadPts,closeLimit1,host);
            fprintf(['%d loops of %d done (%d%%). Avg loop time: %s. ' ...
                     'Expected finish time: %s.\n'],...
                    done,total,floor(done/total*100),...
                    datestr((t-startTime)/(done-preDone),13),datestr(etf));
        end
    end
    results.done=done;
    if doSave
        fprintf('Saving...')
        save(fileName,'results')
        fprintf('done.\n');
    end
end
