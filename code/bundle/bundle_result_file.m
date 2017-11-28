function COP=bundle_result_file(s,e,f,COP)
%BUNDLE_RESULT_FILE Generate result file of bundle run.
%
%   BUNDLE_RESULT_FILE(S,E,F), where S and E are BUNDLE return files and
%   F is a string, writes a text result file to the file F. The text
%   result file contain information about the project, the status of the
%   estimation process, and the quality of the result.
%
%   BUNDLE_RESULT_FILE(S,E,F,COP) supplies a pre-computed OP covariance
%   matrix COP.
%
%   COP=... returns the computed OP covariance matrix.
%
%See also: BUNDLE, BUNDLE_COV. 


fid=fopen(f,'wt');
if fid<0
    error('DBAT:bundle_result_file:fileError',...
          ['Failed to open file ''',f,''' for writing.']);
end

% Signal high correlation above this value.
corrThreshold=0.95;
% Signal low significance below this value.
sigThreshold=0.95;

% Header info
fprintf(fid,'Damped Bundle Adjustment Toolbox result file\n');
p=repmat(' ',1,3);
fprintf(fid,[p,'Project Name: %s\n'],s.title);

fprintf(fid,[p,'Problems and suggestions:\n']);
fprintf(fid,[p,p,'Project Problems: Not evaluated\n']);

% Check if we have any large correlations.
[iio,jio,kio,vio,CIO]=high_io_correlations(s,e,corrThreshold);
[ieo,jeo,keo,veo,CEO]=high_eo_correlations(s,e,corrThreshold);
if nargin<4
    COP=bundle_cov(s,e,'COP');
end
OPstd=sqrt(reshape(full(diag(COP)),3,[]));

[iop,jop,kop,vop]=high_op_correlations(s,e,corrThreshold,COP);
% Compute p values for distortion parameters.
[pk,pp]=test_distortion_params(s,e);
n=any(iio)+any(ieo)+any(iop)+any([pk;pp]<sigThreshold);

fprintf(fid,[p,p,'Problems related to the processing: (%d)\n'],n);

if any(iio)
    fprintf(fid,[p,p,p,'One or more of the camera parameter ' ...
                 'has a high correlation (see below).\n']);
end
if any(ieo)
    fprintf(fid,[p,p,p,'One or more of the camera station parameters ' ...
                 'has a high correlation (see below).\n']);
end
if any(iop)
    fprintf(fid,[p,p,p,'One or more of the object point coordinates ' ...
                 'has a high correlation.\n']);
end
if any([pk;pp]<sigThreshold)
    fprintf(fid,[p,p,p,'One or more estimated lens distortion coefficients ' ...
                 'failed significance test (see below).\n']);
end    

% Info about last bundle.
fprintf(fid,[p,'Information from last bundle\n']);

if e.code==0
    status='OK';
else
    status='fail';
end
hostname=getenv('HOST');
if isempty(hostname)
    hostname='<unknown>';
end
% Values are {name,format string,value}
values={'Last Bundle Run:','%s',e.dateStamp,
        'DBAT version:','%s',e.version,
        'MATLAB version:','%s',version,
        'Host system:','%s',computer,
        'Host name:','%s',hostname,
        'Status:','%s (%d)',{status,e.code},
        'Sigma0:','%g',e.s0,
        'Sigma0 (pixels):','%g',e.s0*s.prior.sigmas(1)};
pretty_print(fid,repmat(p,1,2),values);
    
fprintf(fid,[p,p,'Processing options:\n']);

offon={'off','on'};

values={'Orientation:','%s','on',
        'Global optimization:','%s','on',
        'Calibration:','%s',offon{double(any(s.estIO(:)))+1},
        'Constraints:','%s','off',
        'Maximum # of iterations:','%d',e.maxIter,
        'Convergence tolerance:','%g',e.convTol,
        'Singular test:','%s',offon{double(e.singularTest)+1},
        'Chirality veto:','%s',offon{double(e.chirality)+1},
        'Damping:','%s',e.damping.name,
        'Camera unit (cu):','%s',s.camUnit,
        'Object space unit (ou):','%s',s.objUnit,
        'Initial value comment:','%s',s.x0desc};

pretty_print(fid,repmat(p,1,3),values);

fprintf(fid,[p,p,'Total error:\n']);

values={
    'Number of stages:','%d',1,
    'Number of iterations:','%d',e.usedIters,
    'First error:','%g',e.res(1),
    'Last error:','%g',e.res(end),
    'Execution time (s):','%g',e.time
    };

pretty_print(fid,repmat(p,1,3),values);

fprintf(fid,[p,p,'Lens distortion models:\n']);
if all(s.IOdistModel==1)
    fprintf(fid,[p,p,p,'Backward (Photogrammetry)\n']);
elseif all(s.IOdistModel==-1)
    fprintf(fid,[p,p,p,'Forward (Computer Vision)\n']);
else
    fprintf(fid,[p,p,p,'Mixed Forward/Backward\n']);
end
    
corrStr=sprintf('Correlations over %g%%:',corrThreshold*100);

fprintf(fid,[p,p,'Cameras:\n']);

% Construct a string that indicates what camera parameters are estimated.
selfCal=any(s.estIO,1);
strs={'Xp','Yp','f','K1','K2','K3','P1','P2'};
if all(selfCal)
    allParamCal=all(s.estIO,2);
    anyParamCal=any(s.estIO,2);
    if all(allParamCal==anyParamCal)
        % Any parameters that is estimated in one camera is
        % estimated in all cameras.
        calParams=sprintf('%s ',strs{allParamCal(1:length(strs))});
        selfCalStr=['yes (',strtrim(calParams),')'];
    else
        selfCalStr='yes (mixed parameters)';
    end
elseif ~any(selfCal)
    selfCalStr='no';
else
    selfCalStr='mixed';
end

fprintf(fid,[p,p,p,'Calibration: %s\n'],selfCalStr);

% IO standard deviation. Correlations were computed far above.
ioSigma=reshape(sqrt(diag(CIO)),size(s.IO,1),[]);

% Headers and values to print.
head={'Focal Length','Xp - principal point x','Yp - principal point y',...
      'Fw - format width','Fh - format height',...
      'K1 - radial distortion 1','K2 - radial distortion 2',...
      'K3 - radial distortion 3',...
      'P1 - decentering distortion 1','P2 - decentering distortion 2',...
      'Iw - image width','Ih - image height',...
      'Xr - X resolution','Yr - Y resolution',...
      'Pw - pixel width','Ph - pixel height'};
names={'f','Xp','Yp','Fw','Fh','K1','K2','K3','P1','P2','Iw','Ih','Xr','Yr',...
       'Pw','Ph'};
% Spell out camera unit.
unit0={'cu','cu','cu','cu','cu','cu^(-3)','cu^(-5)','cu^(-7)', ...
       'cu^(-3)','cu^(-3)','px','px','px/cu','px/cu','cu', 'cu'};
unit=strrep(unit0,'cu',s.camUnit);
% Original-to-presentation order mapping and inverse.
rows=[3,1:2,11:12,4:6,7:8,13:14,15:16,17:18]; % Last two are not in real vector.
irows=sparse(rows,1,1:length(rows));

% Flip signs of rows 3,6-10.
S=diag((-1).^double(ismember(1:length(rows),[3,6:10])));

for i=1:length(selfCal)
    % Create fake IO vector.
    sIO=[s.IO;1./s.IO(end-1:end,:)];
    % Extend ioSigma. Need to be fixed when Fw, Fh estimation is implemented.
    ioSigma(end+2,1)=0;
    vals=S*full(sIO(rows,i));
    % Significance values.
    sig=nan(size(rows));
    sigma=zeros(size(rows));
    if selfCal(i)
        sig(6:8)=pk;
        sig(9:10)=pp;
        % Add symmetric correlations.
        iio0=iio;
        iio=[iio;jio];
        jio=[jio;iio0];
        kio=repmat(kio,2,1);
        vio=repmat(vio,2,1);

    
        padLength=length('Significance:');

        sigma=full(ioSigma(rows,i));
    else
        padLength=length('Value:');
    end
    fprintf(fid,[p,p,p,'Camera%d\n'],i);
    fprintf(fid,[p,p,p,p,'Lens distortion model:\n']);
    switch s.IOdistModel(i)
      case 1
        fprintf(fid,[p,p,p,p,p,'Backward (Photogrammetry)\n']);
      case -1
        fprintf(fid,[p,p,p,p,p,'Forward (Computer Vision)\n']);
      otherwise
        fprintf(fid,[p,p,p,p,p,'Unknown'\n]);
    end
    for j=1:length(head)
        fprintf(fid,[p,p,p,p,'%s:\n'],head{j});
        values={'Value:','%g %s',{vals(j),unit{j}}};
        if sigma(j)~=0
            values(end+1,:)={'Deviation:','%.3g %s',{sigma(j),unit{j}}};
        end
        if ~isnan(sig(j))
            values(end+1,:)={'Significance:','p=%.2f',sig(j)};
        end
        highCorr=find(kio==i & iio==rows(j));
        if any(highCorr)
            otherParam=irows(jio(highCorr));
            ss='';
            for kk=1:length(otherParam)
                ss=[ss,sprintf(' %s:%.1f%%,',names{otherParam(kk)},...
                               vio(highCorr(kk))*100)];
            end
            ss(end)='.';
            values(end+1,:)={corrStr,'%s',ss};
        end
        pretty_print(fid,[repmat(p,1,5)],values,padLength,padLength);
    end
    % Print field of view.
    aov=2*atan([s.IO(11:12);norm(s.IO(11:12))]/(2*s.IO(3)))*180/pi;
    fprintf(fid,[p,p,p,'Rated angle of view (h,v,d): (%.0f, %.0f, %.0f) deg\n'],aov);
    % Compute distortion at the sensor corners.
    xx=[1,s.IO(end-3)]+0.5*[-1,1];
    yy=[1,s.IO(end-2)]+0.5*[-1,1];
    corners=[xx([1,1,2,2]);yy([1,2,2,1])];
    xr=corners(1,:)/s.IO(end-1)-s.IO(1);
    yr=corners(2,:)/s.IO(end)+s.IO(2);
    r2=xr.^2+yr.^2;
    xcorrR=xr.*(s.IO(4)*r2+s.IO(5)*r2.^2+s.IO(6)*r2.^3);
    ycorrR=yr.*(s.IO(4)*r2+s.IO(5)*r2.^2+s.IO(6)*r2.^3);
    xcorrT=s.IO(7)*(r2+2*xr.^2)+2*s.IO(8)*xr.*yr;
    ycorrT=s.IO(8)*(r2+2*yr.^2)+2*s.IO(7)*xr.*yr;
    xCorr=xcorrR+xcorrT;
    yCorr=ycorrR+ycorrT;
    mx=max(sqrt(xCorr.^2)+sqrt(yCorr.^2));
    % Length of sensor half-diagonal.
    d=sqrt(s.IO(end-5)^2+s.IO(end-4)^2)/2;
    fprintf(fid,[p,p,p,'Largest distortion: %.2g %s (%.1f px, %.1f%% of half-diagonal)\n'],...
                 mx,s.camUnit,mx*mean(s.IO(end-1:end)),mx/d*100);
end

fprintf(fid,[p,p,'Precisions / Standard Deviations:\n']);

fprintf(fid,[p,p,p,'Photograph Standard Deviations:\n']);

% Get camera station covariances.
% CEO=bundle_cov(s,e,'CEOF');
% Compute corresponding correlations and standard deviations.
[~,eoSigma]=corrmat(CEO,true);
eoSigma=reshape(eoSigma,6,[]);

% Headers and values to print.
head={'Omega','Phi','Kappa','Xc','Yc','Zc'};
unit={'deg','deg','deg','ou','ou','ou'};
names=head;

% Original-to-presentation order mapping and inverse.
rows=[4:6,1:3];
irows=sparse(rows,1,1:length(rows));

% Scaling matrix to degrees.
S=diag([180/pi*ones(1,3),ones(1,3)]);

% Add symmetric correlations.
ieo0=ieo;
ieo=[ieo;jeo];
jeo=[jeo;ieo0];
keo=repmat(keo,2,1);
veo=repmat(veo,2,1);

padLength=length('Deviation:');
for i=1:size(s.EO,2)
    fprintf(fid,[p,p,p,p,'Photo %d: %s\n'],i,s.imNames{i});
    vals=S*s.EO(rows,i);
    sigma=S*eoSigma(rows,i);
    for j=1:6
        fprintf(fid,[p,p,p,p,p,'%s:\n'],head{j});
        values={'Value:','%.6f %s',{vals(j),unit{j}}};
        if sigma(j)~=0
            values(end+1,:)={'Deviation:','%.3g %s',{sigma(j),unit{j}}};
        end
        highCorr=find(keo==i & ieo==rows(j));
        if any(highCorr)
            otherParam=irows(jeo(highCorr));
            ss='';
            for kk=1:length(otherParam)
                ss=[ss,sprintf(' %s:%.1f%%,',names{otherParam(kk)},...
                               veo(highCorr(kk))*100)];
            end
            ss(end)='.';
            values(end+1,:)={corrStr,'%s',ss};
        end
        pretty_print(fid,repmat(p,1,6),values,padLength,padLength);
    end
end

fprintf(fid,[p,'Quality\n']);

fprintf(fid,[p,p,'Photographs\n']);
values={
    'Total number:','%d',length(s.imNames),
    'Numbers used:','%d',size(s.EO,2)
    };
pretty_print(fid,repmat(p,1,3),values);

fprintf(fid,[p,p,'Cameras\n']);
fprintf(fid,[p,p,p,'Total number: %d\n'],size(s.IO,2));
for i=1:size(s.IO,2)
    fprintf(fid,[p,p,p,'Camera%d:\n'],i);

    calStrs={'<not available>','yes'};
    values={
        'Calibration:','%s',calStrs{any(s.estIO(:,i))+1},
        'Number of photos using camera:','%d',nnz(s.cams==i)
        };
    pretty_print(fid,repmat(p,1,4),values);

    % Compute individual and union coverage.
    [c,cr,crr]=coverage(s,find(s.cams==i));
    [uc,ucr,ucrr]=coverage(s,find(s.cams==i),true);
    fprintf(fid,[p,p,p,p,'Photo point coverage:\n']);
    values={
        'Rectangular:','%s',...
        sprintf('%d%%-%d%% (%d%% average, %d%% union)',...
                round(min(cr*100)),round(max(cr*100)),round(mean(cr*100)),...
                round(ucr*100)),
        'Convex hull:','%s',...
        sprintf('%d%%-%d%% (%d%% average, %d%% union)',...
                round(min(c*100)),round(max(c*100)),round(mean(c*100)),...
                round(uc*100)),
        'Radial:','%s',...
        sprintf('%d%%-%d%% (%d%% average, %d%% union)',...
                round(min(crr*100)),round(max(crr*100)),round(mean(crr*100)),...
                round(ucrr*100));
        };
    pretty_print(fid,repmat(p,1,5),values);
end

fprintf(fid,[p,p,'Photo Coverage\n']);
fprintf(fid,[p,p,p,'References points outside calibrated region:\n']);
if any(s.estIO(:))
    fprintf(fid,[p,p,p,p,'none\n']);
else
    fprintf(fid,[p,p,p,p,'<not available>\n']);
end

fprintf(fid,[p,p,'Point Measurements\n']);
fprintf(fid,[p,p,p,'Number of control pts: %d\n'],nnz(s.isCtrl));
fprintf(fid,[p,p,p,'Number of check pts: %d\n'],nnz(s.isCheck));
fprintf(fid,[p,p,p,'Number of object pts: %d\n'],nnz(~s.isCtrl & ~s.isCheck));

if any(s.isCtrl)
    CPrays=full(sum(s.vis(s.isCtrl,:),2));
    n0=nnz(CPrays==0);
    CPrays=CPrays(CPrays~=0);
    mn=min(CPrays);
    mx=max(CPrays);
    avg=mean(CPrays);
    if n0>0
        fprintf(fid,[p,p,p,'CP ray count: %dx0, %d-%d (%.1f avg)\n'],n0,mn,mx,avg);
    else
        fprintf(fid,[p,p,p,'CP ray count: %d-%d (%.1f avg)\n'],mn,mx,avg);
    end
    
    % Ray count histogram.
    cpRayHist=ihist(sum(s.vis(s.isCtrl,:),2)+1);
    i=find(cpRayHist);
    for j=1:length(i)
        fprintf(fid,[p,p,p,p,'%d points with %d rays.\n'],...
                full(cpRayHist(i(j))),i(j)-1);
    end
else
    fprintf(fid,[p,p,p,'CP ray count: -\n']);
end

if any(s.isCheck)
    CCPrays=full(sum(s.vis(s.isCheck,:),2));
    mn=min(CCPrays);
    mx=max(CCPrays);
    avg=mean(CCPrays);
    fprintf(fid,[p,p,p,'CCP ray count: %d-%d (%.1f avg)\n'],mn,mx,avg);

    % Ray count histogram.
    cpRayHist=ihist(sum(s.vis(s.isCheck,:),2)+1);
    i=find(cpRayHist);
    for j=1:length(i)
        fprintf(fid,[p,p,p,p,'%d points with %d rays.\n'],...
                full(cpRayHist(i(j))),i(j)-1);
    end
else
    fprintf(fid,[p,p,p,'CCP ray count: -\n']);
end

if any(~s.isCtrl)
    OPrays=full(sum(s.vis(~s.isCtrl,:),2));
    mn=min(OPrays);
    mx=max(OPrays);
    avg=mean(OPrays);
    fprintf(fid,[p,p,p,'OP ray count: %d-%d (%.1f avg)\n'],mn,mx,avg);
    % Ray count histogram.
    opRayHist=ihist(sum(s.vis(~s.isCtrl,:),2)+1);
    i=find(opRayHist);
    for j=1:length(i)
        fprintf(fid,[p,p,p,p,'%d points with %d rays.\n'],...
                full(opRayHist(i(j))),i(j)-1);
    end
else
    fprintf(fid,[p,p,p,'OP ray count: -\n']);
end

% Get all residuals.
[rms,res]=bundle_residuals(s,e);

fprintf(fid,[p,p,'Point Marking Residuals\n']);
fprintf(fid,[p,p,p,'Overall point RMS: %.3f pixels\n'],rms);

fprintf(fid,[p,p,p,'Mark point residuals:\n']);
fprintf(fid,[p,p,p,p,'Maximum: ']);

[mx,i]=max(res(:));
[mxi,mxj]=ind2sub(size(res),i);
fprintf(fid,'%.3f pixels (OP %d on photo %d)\n',full(mx),s.OPid(mxi),mxj);

% Compute averages per object point.
nOP=full(sum(s.vis,2));
sqSumOP=full(sum(res.^2,2));
meanOP=sqrt(sqSumOP./nOP);

fprintf(fid,[p,p,p,'Object point residuals (RMS over all images of a point):\n']);

[mn,mni]=min(meanOP);
[mx,mxi]=max(meanOP);
fprintf(fid,[p,p,p,p,'Minimum: ']);
fprintf(fid,'%.3f pixels (OP %d over %d images)\n',mn,s.OPid(mni),nOP(mni));
fprintf(fid,[p,p,p,p,'Maximum: ']);
fprintf(fid,'%.3f pixels (OP %d over %d images)\n',mx,s.OPid(mxi),nOP(mxi));

% Compute averages per photo.
nPhoto=full(sum(s.vis,1));
sqSumPhoto=full(sum(res.^2,1));
meanPhoto=sqrt(sqSumPhoto./nPhoto);

fprintf(fid,[p,p,p,'Photo residuals (RMS over all points in an image):\n']);

[mn,mni]=min(meanPhoto);
[mx,mxi]=max(meanPhoto);
fprintf(fid,[p,p,p,p,'Minimum: ']);
fprintf(fid,'%.3f pixels (photo %d over %d points)\n',mn,mni,nPhoto(mni));
fprintf(fid,[p,p,p,p,'Maximum: ']);
fprintf(fid,'%.3f pixels (photo %d over %d points)\n',mx,mxi,nPhoto(mxi));

fprintf(fid,[p,p,'Point Precision\n']);

% Variance for each OP.
v=reshape(full(diag(COP)),3,[]);
% Mask fixed OPs.
v(~s.estOP)=nan;
% Total variance.
tVar=sum(v,1);
% Total standard deviation.
tStd=sqrt(tVar);

fprintf(fid,[p,p,p,'Total standard deviation (RMS of X/Y/Z std):\n']);

[mn,mni]=min(tStd);
[mx,mxi]=max(tStd);
fprintf(fid,[p,p,p,p,'Minimum: ']);
fprintf(fid,'%.2g (OP %d)\n',mn,s.OPid(mni));
fprintf(fid,[p,p,p,p,'Maximum: ']);
fprintf(fid,'%.2g (OP %d)\n',mx,s.OPid(mxi));

% Component-wise.
[mx,mxi]=max(sqrt(v),[],2);
for i=1:3
    fprintf(fid,[p,p,p,'Maximum %c standard deviation: %.2g (OP %d)\n'],...
                 abs('X')+i-1,mx(i),s.OPid(mxi(i)));
end

% Points with highest correlations.
fprintf(fid,[p,p,p,'Points with high correlations\n']);
opCorr95=abs(vop)>0.95;
opCorr99=abs(vop)>0.99;
fprintf(fid,[p,p,p,p,'Points with correlation above 95%%: %d\n'],nnz(opCorr95));
fprintf(fid,[p,p,p,p,'Points with correlation above 99%%: %d\n'],nnz(opCorr99));
if nnz(opCorr95)
    fprintf(fid,[p,p,p,p,'Points with highest correlations:\n']);
    [~,i]=sort(abs(vop),'descend');
    printed=[];
    j=1;
    while length(printed)<5 && j<=length(i)
        % Same OP may appear multiple times due to e.g. high X-Z,
        % Y-Z corr. Only print each OP once.
        if ~ismember(kop(i(j)),printed)
            printed(end+1)=kop(i(j));
            fprintf(fid,[p,p,p,p,p,'Points %d: %.2f\n'], kop(i(j)),...
                    100*vop(i(j)));
        end
        j=j+1;
    end
end

fprintf(fid,[p,p,'Point Angles\n']);

% Maximum angle between rays.
a=angles(s,'Computing angles')*180/pi;

fprintf(fid,[p,p,p,'CP\n']);

% Ctrl pts with >0 rays
cpIx=find(s.isCtrl);
cpRays=sum(s.vis(cpIx,:),2);
if any(cpRays==0)
    fprintf(fid,[p,p,p,p,'Ignoring %d CP with 0 rays.\n'],nnz(cpRays==0));
end
% If we have any ctrl pts with rays.
if any(cpIx(cpRays>0))
    aCP=a(s.isCtrl);
    idCP=s.OPid(s.isCtrl);
    labelCP=s.OPlabels(s.isCtrl);
    [mn,mni]=min(aCP);
    [mx,mxi]=max(aCP);
    if isempty(labelCP{mni})
        lmn='';
    else
        lmn=sprintf(', label %s',labelCP{mni});
    end
    if isempty(labelCP{mxi})
        lmx='';
    else
        lmx=sprintf(', label %s',labelCP{mxi});
    end
    fprintf(fid,[p,p,p,p,'Minimum: %.1f degrees (CP %d%s)\n'],mn,idCP(mni),lmn);
    fprintf(fid,[p,p,p,p,'Maximum: %.1f degrees (CP %d%s)\n'],mx,idCP(mxi),lmx);
    fprintf(fid,[p,p,p,p,'Average: %.1f degrees\n'],nanmean(aCP));
else
    fprintf(fid,[p,p,p,p,'Minimum: -\n']);
    fprintf(fid,[p,p,p,p,'Maximum: -\n']);
    fprintf(fid,[p,p,p,p,'Average: -\n']);
end

fprintf(fid,[p,p,p,'CCP\n']);
if any(s.isCheck)
    aCCP=a(s.isCheck);
    idCCP=s.OPid(s.isCheck);
    labelCCP=s.OPlabels(s.isCheck);
    [mn,mni]=min(aCCP);
    [mx,mxi]=max(aCCP);
    if isempty(labelCCP{mni})
        lmn='';
    else
        lmn=sprintf(', label %s',labelCCP{mni});
    end
    if isempty(labelCCP{mxi})
        lmx='';
    else
        lmx=sprintf(', label %s',labelCCP{mxi});
    end
    fprintf(fid,[p,p,p,p,'Minimum: %.1f degrees (CCP %d%s)\n'],mn,idCCP(mni),lmn);
    fprintf(fid,[p,p,p,p,'Maximum: %.1f degrees (CCP %d%s)\n'],mx,idCCP(mxi),lmx);
    fprintf(fid,[p,p,p,p,'Average: %.1f degrees\n'],mean(aCCP));
else
    fprintf(fid,[p,p,p,p,'Minimum: -\n']);
    fprintf(fid,[p,p,p,p,'Maximum: -\n']);
    fprintf(fid,[p,p,p,p,'Average: -\n']);
end

fprintf(fid,[p,p,p,'OP\n']);
isOP=~s.isCtrl & ~s.isCheck;
if any(isOP)
    aOP=a(isOP);
    idOP=s.OPid(isOP);
    [mn,mni]=min(aOP);
    [mx,mxi]=max(aOP);
    fprintf(fid,[p,p,p,p,'Minimum: %.1f degrees (OP %d)\n'],mn,idOP(mni));
    fprintf(fid,[p,p,p,p,'Maximum: %.1f degrees (OP %d)\n'],mx,idOP(mxi));
    fprintf(fid,[p,p,p,p,'Average: %.1f degrees\n'],mean(aOP));
    fprintf(fid,[p,p,p,p,'Smallest angles (ID, angle [deg], vis in ' ...
                 'cameras)\n']);
    [ang,i]=sort(aOP);
    % Get all points with at least third min angle + 10% + 0.1 deg...
    angLimit=ang(min(3,end))*1.1+0.1;
    % ...but no point above 80 degress.
    angLimit=min(angLimit,80);
    % At least 3 points but obviously not more than all.
    nPts=min(max(nnz(ang<angLimit),3),length(ang));
    for j=1:nPts
        camVis=find(s.vis(i(j),:));
        str=sprintf('%4d ',camVis);
        fprintf(fid,[p,p,p,p,p,'%6d: %5.2f (%s)\n'],idOP(i(j)),ang(j),str(1:end-1));
    end
else
    fprintf(fid,[p,p,p,p,'Minimum: -\n']);
    fprintf(fid,[p,p,p,p,'Maximum: -\n']);
    fprintf(fid,[p,p,p,p,'Average: -\n']);
end

fprintf(fid,[p,p,'Ctrl measurements\n']);
if any(s.isCtrl)
    cIx=find(s.isCtrl);
    CPid=s.OPid(cIx);

    fprintf(fid,[p,p,p,'Prior\n']);
    fprintf(fid,[p,p,p,'%6s, %8s, %8s, %8s, %8s, %8s, %8s, %s\n'],...
            'id','x','y','z','stdx','stdy','stdz','label');
    pos0=s.prior.OP(:,cIx);
    std0=s.prior.OPstd(:,cIx);
    for i=1:length(cIx)
        fprintf(fid,[p,p,p,'%6d, %8.3f, %8.3f, %8.3f, %8.3g, ' ...
                     '%8.3g, %8.3g, %s\n'],CPid(i),pos0(:,i),std0(:,i),...
                s.OPlabels{cIx(i)});
    end
    fprintf(fid,[p,p,p,'Posterior\n']);
    fprintf(fid,[p,p,p,'%6s, %8s, %8s, %8s, %8s, %8s, %8s, %4s, %s\n'],...
            'id','x','y','z','stdx','stdy','stdz','rays','label');
    pos1=s.OP(:,cIx);
    std1=OPstd(:,cIx);
    for i=1:length(cIx)
        fprintf(fid,[p,p,p,'%6d, %8.3f, %8.3f, %8.3f, %8.3g, ' ...
                     '%8.3g, %8.3g, %4d, %s\n'],CPid(i),pos1(:,i),std1(:,i),...
                full(sum(s.vis(cIx(i),:))),s.OPlabels{cIx(i)});
    end
    fprintf(fid,[p,p,p,'Diff (pos=abs diff, std=rel diff)\n']);
    fprintf(fid,[p,p,p,'%6s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %4s, %s\n'],...
            'id','x','y','z','xy','xyz','stdx','stdy','stdz','rays','label');
    posd=pos1-pos0;
    stdd=((std1+eps)./(std0+eps)-1)*100;
    for i=1:length(cIx)
        fprintf(fid,[p,p,p,'%6d, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %7.1f%%, ' ...
        '%7.1f%%, %7.1f%%, %4d, %s\n'],CPid(i),posd(:,i),norm(posd(1:2,i)),...
                norm(posd(:,i)),stdd(:,i),...
                full(sum(s.vis(cIx(i),:))),s.OPlabels{cIx(i)});
    end

    fprintf(fid,[p,p,p,'Ctrl point delta\n']);
    diffNorm=sqrt(sum(posd.^2));
    [mx,i]=max(diffNorm);
    maxId=CPid(i);
    maxLabel=s.OPlabels{cIx(i)};
    fprintf(fid,[p,p,p,p,'Max: %.3f ou (%s, pt %d)\n'],mx,maxLabel,maxId);
    fprintf(fid,[p,p,p,p,'Max X,Y,Z\n']);
    for i=1:3
        [mx,j]=max(abs(posd(i,:)));
        fprintf(fid,[p,p,p,p,p,'%c: %.3f ou (%s, pt %d)\n'],'X'-1+i,mx,...
                s.OPlabels{cIx(j)},CPid(j));
    end
    fprintf(fid,[p,p,p,p,'RMS: %.3f ou (from %d items)\n'],...
            sqrt(mean(diffNorm.^2)),length(diffNorm));
end

fprintf(fid,[p,p,'Check measurements\n']);
if any(s.isCheck)
    cIx=find(s.isCheck);
    CPid=s.OPid(cIx);

    fprintf(fid,[p,p,p,'Prior\n']);
    fprintf(fid,[p,p,p,'%6s, %8s, %8s, %8s, %8s, %8s, %8s, %s\n'],...
            'id','x','y','z','stdx','stdy','stdz','label');
    pos0=s.prior.CCP(:,cIx);
    std0=s.prior.CCPstd(:,cIx);
    for i=1:length(cIx)
        fprintf(fid,[p,p,p,'%6d, %8.3f, %8.3f, %8.3f, %8.3g, ' ...
                     '%8.3g, %8.3g, %s\n'],CPid(i),pos0(:,i),std0(:,i),...
                s.OPlabels{cIx(i)});
    end
    fprintf(fid,[p,p,p,'Posterior\n']);
    fprintf(fid,[p,p,p,'%6s, %8s, %8s, %8s, %8s, %8s, %8s, %4s, %s\n'],...
            'id','x','y','z','stdx','stdy','stdz','rays','label');
    pos1=s.OP(:,cIx);
    std1=OPstd(:,cIx);
    for i=1:length(cIx)
        fprintf(fid,[p,p,p,'%6d, %8.3f, %8.3f, %8.3f, %8.3g, ' ...
                     '%8.3g, %8.3g, %4d, %s\n'],CPid(i),pos1(:,i),std1(:,i),...
                full(sum(s.vis(cIx(i),:))),s.OPlabels{cIx(i)});
    end
    fprintf(fid,[p,p,p,'Diff (pos=abs diff, std=rel diff)\n']);
    fprintf(fid,[p,p,p,'%6s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %4s, %s\n'],...
            'id','x','y','z','xy','xyz','stdx','stdy','stdz','rays','label');
    posd=pos1-pos0;
    stdd=((std1+eps)./(std0+eps)-1)*100;
    for i=1:length(cIx)
        fprintf(fid,[p,p,p,'%6d, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %7.1f%%, ' ...
        '%7.1f%%, %7.1f%%, %4d, %s\n'],CPid(i),posd(:,i),norm(posd(1:2,i)),...
                norm(posd(:,i)),stdd(:,i),...
                full(sum(s.vis(cIx(i),:))),s.OPlabels{cIx(i)});
    end

    fprintf(fid,[p,p,p,'Check point delta\n']);
    diffNorm=sqrt(sum(posd.^2));
    [mx,i]=max(diffNorm);
    maxId=CPid(i);
    maxLabel=s.OPlabels{cIx(i)};
    fprintf(fid,[p,p,p,p,'Max: %.3f ou (%s, pt %d)\n'],mx,maxLabel,maxId);
    fprintf(fid,[p,p,p,p,'Max X,Y,Z\n']);
    for i=1:3
        [mx,j]=max(abs(posd(i,:)));
        fprintf(fid,[p,p,p,p,p,'%c: %.3f ou (%s, pt %d)\n'],'X'-1+i,mx,...
                s.OPlabels{cIx(j)},CPid(j));
    end
    fprintf(fid,[p,p,p,p,'RMS: %.3f ou (from %d items)\n'],...
            sqrt(mean(diffNorm.^2)),length(diffNorm));
end

fclose(fid);


function pretty_print(fid,prefix,values,minLen,maxLen)
% Pretty-prints the values in values with prefix.

if nargin<4, minLen=inf; end
if nargin<5, maxLen=-inf; end

% Determine padding.
nameLen=cellfun(@length,values(:,1));
pad=max(max(min(minLen,max(nameLen)),maxLen)+1-nameLen,0);
for i=1:size(values,1)
    val=values{i,3};
    if ~iscell(val), val={val}; end
    fprintf(fid,[prefix,'%s%s',values{i,2},'\n'],values{i,1},...
            repmat(' ',1,pad(i)),val{:});
end
