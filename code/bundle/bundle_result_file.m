function [COP,s]=bundle_result_file(s,e,f,COP)
%BUNDLE_RESULT_FILE Generate result file of bundle run.
%
%   BUNDLE_RESULT_FILE(S,E,F), where S and E are BUNDLE return structs
%   and F is a string, writes a text result file to the file F. The
%   text result file contain information about the project, the status
%   of the estimation process, and the quality of the result.
%
%   BUNDLE_RESULT_FILE(S,E,F,COP) supplies a pre-computed OP covariance
%   matrix COP.
%
%   [COP,S]=... returns the computed OP covariance matrix and the
%   updated DBAT structure.
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
p2=repmat(p,1,2);
p3=repmat(p,1,3);
p4=repmat(p,1,4);
p5=repmat(p,1,5);
p6=repmat(p,1,6);

fprintf(fid,[p,'Project\n']);

fprintf(fid,[p2,'Name        : %s\n'],s.proj.title);

if ~isempty(s.proj.fileName)
    str=s.proj.fileName;
    str=strrep(str,[dbatroot,filesep],['$DBATROOT',filesep]);
    fprintf(fid,[p2,'File name   : %s\n'],str);
end

if ~isempty(s.proj.cptFile)
    str=s.proj.cptFile;
    str=strrep(str,[dbatroot,filesep],['$DBATROOT',filesep]);
    fprintf(fid,[p2,'Ctrl pt file: %s\n'],str);
end

if ~isempty(s.proj.EOfile)
    str=s.proj.EOfile;
    str=strrep(str,[dbatroot,filesep],['$DBATROOT',filesep]);
    fprintf(fid,[p2,'EO file     : %s\n'],str);
end
    
    
fprintf(fid,[p,'Problems and suggestions:\n']);
fprintf(fid,[p2,'Project Problems:\n']);
if isempty(e.weakness.structural)
    fprintf(fid,[p3,'Structural rank: ok.\n']);
else
    fprintf(fid,[p3,'Structural rank: %d (deficiency: %d)\n'],...
            e.weakness.structural.rank,e.weakness.structural.deficiency);
    fprintf(fid,[p4,'DMPERM suggests the ' ...
                        'following parameters have problems:\n']);
    for i=1:length(e.weakness.structural.suspectedParams)
        fprintf(fid,[p5,'%s\n'],e.weakness.structural.suspectedParams{i});
    end
end
if isempty(e.weakness.numerical) || e.weakness.numerical.deficiency==0
    fprintf(fid,[p3,'Numerical rank: ok.\n']);
elseif isnan(e.weakness.numerical.rank)
    fprintf(fid,[p3,'Numerical rank: not tested.\n']);
else
    fprintf(fid,[p3,'Numerical rank: %d (deficiency: %d)\n'],...
            e.weakness.numerical.rank,e.weakness.numerical.deficiency);
    fprintf(fid,[p4,'Null-space suggest the following parameters ' ...
                 'are part of the problem:\n']);
    
    for i=1:length(e.weakness.numerical.suspectedParams)
        fprintf(fid,[p5,'Vector %d (eigenvalue %g):\n'],...
                i,e.weakness.numerical.d(i));
        sp=e.weakness.numerical.suspectedParams{i};
        for j=1:length(sp.values)
            fprintf(fid,[p6,'(%s, %.3g)\n'],sp.params{j},sp.values(j));
        end
    end
end

% Check if we have any large correlations.
[iio,jio,vio,CIO]=high_io_correlations(s,e,corrThreshold,true);
[ieo,jeo,keo,veo,CEO]=high_eo_correlations(s,e,corrThreshold);
if nargin<4
    COP=bundle_cov(s,e,'COP');
end
[iop,jop,kop,vop]=high_op_correlations(s,e,corrThreshold,COP);

% Compute posterior standard devaitions.
s.post.std.IO=sqrt(reshape(full(diag(CIO)),size(s.IO.val)));
s.post.std.EO=sqrt(reshape(full(diag(CEO)),size(s.EO.val)));
s.post.std.OP=sqrt(reshape(full(diag(COP)),size(s.OP.val)));

% Compute p values for distortion parameters.
[pk,pp,pb,pkc]=test_distortion_params(s,e);
n=double(any(vio))+double(any(veo))+double(any(vop))+ ...
  double(any(any([pk;pp;pb]<sigThreshold)))+double(e.code~=0);

fprintf(fid,[p2,'Problems related to the processing: (%d)\n'],n);

if e.code~=0
    fprintf(fid,[p3,'Bundle failed with code %d (see below for details).\n'],...
            e.code);
end
if any(iio)
    fprintf(fid,[p3,'One or more of the camera parameter ' ...
                 'has a high correlation (see below).\n']);
end
if any(ieo)
    fprintf(fid,[p3,'One or more of the camera station parameters ' ...
                 'has a high correlation (see below).\n']);
end
if any(iop)
    fprintf(fid,[p3,'One or more of the object point coordinates ' ...
                 'has a high correlation.\n']);
end
if any(any([pk;pp;pb]<sigThreshold))
    fprintf(fid,[p3,'One or more estimated lens and/or affine distortion coefficients ' ...
                 'failed significance test (see below).\n']);
end    

% Info about last bundle.
fprintf(fid,[p,'Information from last bundle\n']);

if e.code==0
    status='OK';
else
    msgs={'Too many iterations','Normal matrix is singular',...
          'No step length found by the line search',...
          'Normal matrix is structurally rank deficient'};
    if abs(e.code)<=length(msgs)
        status=sprintf('fail (code %d: %s)',e.code,msgs{abs(e.code)});
    else
        status=sprintf('fail (code %d: unknown code)',e.code);
    end
end
hostname=getenv('HOST');
if isempty(hostname)
    hostname='<unknown>';
end

numParams=sprintf('%d (%d IO, %d EO, %d OP)',e.numParams,...
                  nnz(s.IO.struct.leading),nnz(s.EO.struct.leading),...
                  nnz(s.bundle.est.OP));
numObs=sprintf('%d (%d IP, %d IO, %d EO, %d OP)',...
               e.numObs,length(s.post.res.ix.IP),...
               length(s.post.res.ix.IO),...
               length(s.post.res.ix.EO),...
               length(s.post.res.ix.OP));

[comp,mxSize,endian]=computer;
compStr=sprintf('%s (endian=%s, max #elems=%d)',comp,endian,mxSize);

% Values are {name,format string,value}
values={'Last Bundle Run:','%s',e.dateStamp,
        'DBAT version:','%s',e.version,
        'MATLAB version:','%s',version,
        'Host system:','%s',compStr,
        'Host name:','%s',hostname,
        'Status:','%s',status,
        'Sigma0:','%g',e.s0,
        'Sigma0 (pixels):','%g',s.post.sigmas(1),
        'Redundancy','%d',e.redundancy,
        'Number of params:','%s',numParams,
        'Number of observations:','%s',numObs};
pretty_print(fid,p2,values);
    
fprintf(fid,[p2,'Processing options:\n']);

offon={'off','on'};
relabs={'relative','absolute'};

values={'Orientation:','%s','on',
        'Global optimization:','%s','on',
        'Calibration:','%s',offon{double(any(s.bundle.est.IO(:)))+1},
        'Constraints:','%s','off',
        'Maximum # of iterations:','%d',e.maxIter,
        'Convergence tolerance:','%g',e.convTol,
        'Termination criteria:','%s',relabs{double(e.absTerm)+1},
        'Singular test:','%s',offon{double(e.singularTest)+1},
        'Chirality veto:','%s',offon{double(e.chirality)+1},
        'Damping:','%s',e.damping.name,
        'Camera unit (cu):','%s',s.IO.model.camUnit,
        'Object space unit (ou):','%s',s.proj.objUnit,
        'Initial value comment:','%s',s.proj.x0desc};

pretty_print(fid,p3,values);

fprintf(fid,[p2,'Total error:\n']);

values={
    'Number of stages:','%d',1,
    'Number of iterations:','%d',e.usedIters,
    'First error:','%g',e.res(1),
    'Last error:','%g',e.res(end),
    'Execution time (s):','%g',e.time
    };

pretty_print(fid,p3,values);

fprintf(fid,[p2,'Lens distortion models:\n']);
distModel=unique(s.IO.model.distModel);
if isscalar(distModel) && distModel>0
    fprintf(fid,[p3,'Backward (Photogrammetry) model %d\n'],distModel);
elseif isscalar(distModel) && distModel<0
    fprintf(fid,[p3,'Forward (Computer Vision) model %d\n'],-distModel);
else
    fprintf(fid,[p3,'Mixed Forward/Backward\n']);
end
    
corrStr=sprintf('Correlations over %g%%:',corrThreshold*100);

fprintf(fid,[p2,'Cameras:\n']);

% Construct a string that indicates what camera parameters are estimated.
selfCal=any(s.bundle.est.IO,1);
strs=s.IO.type(:,1);
if all(selfCal)
    allParamCal=all(s.bundle.est.IO,2);
    anyParamCal=any(s.bundle.est.IO,2);
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

fprintf(fid,[p3,'Calibration: %s\n'],selfCalStr);

% Headers and values to print.
names={'cc','px','py','fw','fh','K1','K2','K3','P1','P2','as', ...
       'sk','iw','ih','xr','yr', 'pw','ph'};
head={'Camera Constant','px - principal point x','py - principal point y',...
      'Format width','Format height',...
      'K1 - radial distortion 1','K2 - radial distortion 2',...
      'K3 - radial distortion 3',...
      'P1 - decentering distortion 1','P2 - decentering distortion 2',...
      'as - off-unit aspect parameter','sk - skew',...
      'Image width','Image height',...
      'X resolution','Y resolution',...
      'Pixel width','Pixel height'};

% Spell out camera unit.
unit0={'cu','cu','cu','cu','cu','cu^(-3)','cu^(-5)','cu^(-7)', ...
       'cu^(-3)','cu^(-3)','','','px','px','px/cu','px/cu','cu', 'cu'};
unit=strrep(unit0,'cu',s.IO.model.camUnit);

% Types of data stored in the IOsensorData array.
sensorTypes={'fw','fh','iw','ih','xr','yr','pw','ph'};

% Array with IO data, std and significance values.
IOsensorData=[s.post.sensor.ssSize;
              s.post.sensor.imSize;
              s.post.sensor.imSize./s.post.sensor.ssSize;
              s.post.sensor.pxSize];
IOdata=s.IO.val;
IOtype=s.IO.type;
IOdataStd=s.post.std.IO;
IOdataSig=[nan(3,size(s.IO.val,2));
           pb;
           pk;
           pp];
IOdataCumSig=[nan(3,size(s.IO.val,2));
           nan(size(pb));
           pk;
           nan(size(pp))];

% Original-to-presentation order mapping and inverse.
rows=[1:3,-(1:2),5+(1:s.IO.model.nK+s.IO.model.nP),4:5,-(3:8)];

% Flip signs of py.
IOdata([3,6:end],:)=-IOdata([3,6:end],:);

% For each camera with parameters that were estimated
for i=find(s.IO.struct.uniq)
    if selfCal(i)
        % Correlations [i,k1,j,k2,v]. IO parameter i in camera k1
        % is correlated with IO parameter j in camera k2 with
        % correlation v.
        cc=[iio,jio,vio;
            jio,iio,vio];
        padLength=length('Significance:');
    else
        cc=zeros(0,5);
        padLength=length('Value:');
    end
    if s.IO.struct.isSimple(i)
        fprintf(fid,[p3,'Camera%d (simple)\n'],s.IO.struct.no(i));
    else
        fprintf(fid,[p3,'Camera%d (mixed)\n'],s.IO.struct.no(i));
    end
    fprintf(fid,[p4,'Lens distortion model:\n']);
    if s.IO.model.distModel(i)>0
        fprintf(fid,[p5,'Backward (Photogrammetry) model %d\n'],...
                s.IO.model.distModel(i));
    else
        fprintf(fid,[p5,'Forward (Computer Vision) model %d\n'],...
                -s.IO.model.distModel(i));
    end
    for j=1:length(head)
        jj=rows(j);
        if jj>0
            % Camera data
            val=IOdata(jj,i);
            sigma=IOdataStd(jj,i);
            sig=IOdataSig(jj,i);
            cumSig=IOdataCumSig(jj,i);
        else
            % Sensor data
            val=IOsensorData(abs(jj),i);
            sigma=nan;
            sig=nan;
            cumSig=nan;
        end
        fprintf(fid,[p4,'%s:\n'],head{j});
        values={'Value:','%g %s',{val,unit{j}}};
        if ~isnan(sigma) && sigma~=0
            values(end+1,:)={'Deviation:','%.3g %s',{sigma,unit{j}}}; %#ok<AGROW>
        end
        if ~isnan(sig)
            values(end+1,:)={'Significance:','p=%.2f',sig}; %#ok<AGROW>
        end
        if ~isnan(cumSig)
            values(end+1,:)={'Cumulative significance:','p=%.2f',cumSig}; %#ok<AGROW>
        end
        highCorr=find(cc(:,2)==i & cc(:,1)==jj);
        if any(highCorr)
            otherParam=cc(highCorr,3);
            otherCam=cc(highCorr,4);
            corrVal=cc(highCorr,5);
            ss='';
            for kk=1:length(otherParam)
                if otherCam(kk)==i
                    % Same camera
                    ss=[ss,sprintf(' %s:%.1f%%,',IOtype{otherParam(kk),i},...
                                   corrVal(kk)*100)];
                else
                    % Other camera
                    ss=[ss,sprintf(' %s(cam%d):%.1f%%,',...
                                   IOtype{otherParam(kk),i},...
                                   otherCam(kk),corrVal(kk)*100)];
                end
            end
            ss(end)='.';
            values(end+1,:)={corrStr,'%s',ss}; %#ok<AGROW>
        end
        pretty_print(fid,p5,values,padLength,padLength);
    end
    % Sensor width, height, diagonal
    sensWHD=[s.IO.sensor.ssSize(:,i);norm(s.IO.sensor.ssSize(:,i))];
    % Print field of view.
    aov=2*atan(sensWHD/(2*s.IO.val(1,i)))*180/pi;
    fprintf(fid,[p3,'Rated angle of view (h,v,d): (%.0f, %.0f, %.0f) deg\n'],...
            aov);
    % Compute distortion at the sensor corners.
    xx=[1,s.IO.sensor.imSize(1,i)]+0.5*[-1,1];
    yy=[1,s.IO.sensor.imSize(2,i)]+0.5*[-1,1];
    corners=[xx([1,1,2,2]);yy([1,2,2,1])];
    xr=corners(1,:)*s.IO.sensor.pxSize(1,i)-s.IO.val(2,i);
    yr=corners(2,:)*s.IO.sensor.pxSize(2,i)+s.IO.val(3,i);
    r2=xr.^2+yr.^2;
    xcorrR=xr.*(s.IO.val(6,i)*r2+s.IO.val(7,i)*r2.^2+s.IO.val(8,i)*r2.^3);
    ycorrR=yr.*(s.IO.val(6,i)*r2+s.IO.val(7,i)*r2.^2+s.IO.val(8,i)*r2.^3);
    xcorrT=s.IO.val(9,i)*(r2+2*xr.^2)+2*s.IO.val(9,i)*xr.*yr;
    ycorrT=s.IO.val(10,i)*(r2+2*yr.^2)+2*s.IO.val(10,i)*xr.*yr;
    xCorr=xcorrR+xcorrT;
    yCorr=ycorrR+ycorrT;
    mx=max(sqrt(xCorr.^2)+sqrt(yCorr.^2));
    % Length of sensor half-diagonal.
    d=sensWHD(3)/2;
    fprintf(fid,[p3,'Largest distortion: %.2g %s (%.1f px, %.1f%% of half-diagonal)\n'],...
                 mx,s.IO.model.camUnit,mx/s.IO.sensor.pxSize(1,i),mx/d*100);
end

fprintf(fid,[p2,'Precisions / Standard Deviations:\n']);

fprintf(fid,[p3,'Photograph Standard Deviations:\n']);

% Get camera station covariances.
% CEO=bundle_cov(s,e,'CEOF');
% Compute corresponding correlations and standard deviations.
eoSigma=s.post.std.EO;

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
for i=1:size(s.EO.val,2)
    fprintf(fid,[p4,'Photo %d: %s\n'],i,s.EO.name{i});
    vals=S*s.EO.val(rows,i);
    sigma=S*eoSigma(rows,i);
    for j=1:6
        fprintf(fid,[p5,'%s:\n'],head{j});
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
        pretty_print(fid,p6,values,padLength,padLength);
    end
end

fprintf(fid,[p,'Quality\n']);

fprintf(fid,[p2,'Photographs\n']);
values={
    'Total number:','%d',length(s.EO.name),
    'Numbers used:','%d',size(s.EO.val,2)
    };
pretty_print(fid,p3,values);

fprintf(fid,[p2,'Cameras\n']);
fprintf(fid,[p3,'Total number: %d (%d simple, %d mixed)\n'],...
        nnz(s.IO.struct.uniq),nnz(s.IO.struct.uniq & s.IO.struct.isSimple),...
        nnz(s.IO.struct.uniq & ~s.IO.struct.isSimple));
for i=find(s.IO.struct.uniq)
    fprintf(fid,[p3,'Camera%d:\n'],s.IO.struct.no(i));

    calStrs={'<not available>','yes'};
    values={
        'Calibration:','%s',calStrs{any(s.bundle.est.IO(:,i))+1},
        'Number of photos using camera:','%d',nnz(s.IO.struct.no==s.IO.struct.no(i))
        };
    pretty_print(fid,p4,values);

    % Compute individual and union coverage.
    [c,cr,crr]=coverage(s,find(s.IO.struct.no==s.IO.struct.no(i)));
    [uc,ucr,ucrr]=coverage(s,find(s.IO.struct.no==s.IO.struct.no(i)),true);
    fprintf(fid,[p4,'Photo point coverage:\n']);
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
    pretty_print(fid,p5,values);
end

fprintf(fid,[p2,'Photo Coverage\n']);
fprintf(fid,[p3,'Reference points outside calibrated region:\n']);
for i=find(s.IO.struct.uniq)
    if any(s.bundle.est.IO(:,i))
        fprintf(fid,[p4,'Camera %d: none\n'],s.IO.struct.no(i));
    else
        fprintf(fid,[p4,'Camera %d: <not available>\n'],s.IO.struct.no(i));
    end
end

fprintf(fid,[p2,'Point Measurements\n']);
fprintf(fid,[p3,'Number of control pts: %d\n'],nnz(s.prior.OP.isCtrl));
fprintf(fid,[p3,'Number of check pts: %d\n'],nnz(s.prior.OP.isCheck));
fprintf(fid,[p3,'Number of object pts: %d\n'],nnz(~s.prior.OP.isCtrl & ~s.prior.OP.isCheck));

if any(s.prior.OP.isCtrl)
    CPrays=full(sum(s.IP.vis(s.prior.OP.isCtrl,:),2));
    n0=nnz(CPrays==0);
    CPrays=CPrays(CPrays~=0);
    mn=min(CPrays);
    mx=max(CPrays);
    avg=mean(CPrays);
    if n0>0
        fprintf(fid,[p3,'CP ray count: %dx0, %d-%d (%.1f avg)\n'],n0,mn,mx,avg);
    else
        fprintf(fid,[p3,'CP ray count: %d-%d (%.1f avg)\n'],mn,mx,avg);
    end
    
    % Ray count histogram.
    cpRayHist=ihist(sum(s.IP.vis(s.prior.OP.isCtrl,:),2)+1);
    i=find(cpRayHist);
    for j=1:length(i)
        fprintf(fid,[p4,'%d points with %d rays.\n'],...
                full(cpRayHist(i(j))),i(j)-1);
    end
else
    fprintf(fid,[p3,'CP ray count: -\n']);
end

if any(s.prior.OP.isCheck)
    CCPrays=full(sum(s.IP.vis(s.prior.OP.isCheck,:),2));
    mn=min(CCPrays);
    mx=max(CCPrays);
    avg=mean(CCPrays);
    fprintf(fid,[p3,'CCP ray count: %d-%d (%.1f avg)\n'],mn,mx,avg);

    % Ray count histogram.
    cpRayHist=ihist(sum(s.IP.vis(s.prior.OP.isCheck,:),2)+1);
    i=find(cpRayHist);
    for j=1:length(i)
        fprintf(fid,[p4,'%d points with %d rays.\n'],...
                full(cpRayHist(i(j))),i(j)-1);
    end
else
    fprintf(fid,[p3,'CCP ray count: -\n']);
end

if any(~s.prior.OP.isCtrl)
    OPrays=full(sum(s.IP.vis(~s.prior.OP.isCtrl,:),2));
    mn=min(OPrays);
    mx=max(OPrays);
    avg=mean(OPrays);
    fprintf(fid,[p3,'OP ray count: %d-%d (%.1f avg)\n'],mn,mx,avg);
    % Ray count histogram.
    opRayHist=ihist(sum(s.IP.vis(~s.prior.OP.isCtrl,:),2)+1);
    i=find(opRayHist);
    for j=1:length(i)
        fprintf(fid,[p4,'%d points with %d rays.\n'],...
                full(opRayHist(i(j))),i(j)-1);
    end
else
    fprintf(fid,[p3,'OP ray count: -\n']);
end

% Get all residuals.
[rms,res]=bundle_residuals(s,e);

fprintf(fid,[p2,'Point Marking Residuals\n']);
fprintf(fid,[p3,'Overall point RMS: %.3f pixels\n'],rms);

fprintf(fid,[p3,'Mark point residuals:\n']);
fprintf(fid,[p4,'Maximum: ']);

[mx,i]=max(res(:));
[mxi,mxj]=ind2sub(size(res),i);
fprintf(fid,'%.3f pixels (OP %d on photo %d)\n',full(mx),s.OP.id(mxi),mxj);

% Compute averages per object point.
nOP=full(sum(s.IP.vis,2));
sqSumOP=full(sum(res.^2,2));
meanOP=sqrt(sqSumOP./nOP);

fprintf(fid,[p3,'Object point residuals (RMS over all images of a point):\n']);

[mn,mni]=min(meanOP);
[mx,mxi]=max(meanOP);
fprintf(fid,[p4,'Minimum: ']);
fprintf(fid,'%.3f pixels (OP %d over %d images)\n',mn,s.OP.id(mni),nOP(mni));
fprintf(fid,[p4,'Maximum: ']);
fprintf(fid,'%.3f pixels (OP %d over %d images)\n',mx,s.OP.id(mxi),nOP(mxi));

% Compute averages per photo.
nPhoto=full(sum(s.IP.vis,1));
sqSumPhoto=full(sum(res.^2,1));
meanPhoto=sqrt(sqSumPhoto./nPhoto);

fprintf(fid,[p3,'Photo residuals (RMS over all points in an image):\n']);

[mn,mni]=min(meanPhoto);
[mx,mxi]=max(meanPhoto);
fprintf(fid,[p4,'Minimum: ']);
fprintf(fid,'%.3f pixels (photo %d over %d points)\n',mn,mni,nPhoto(mni));
fprintf(fid,[p4,'Maximum: ']);
fprintf(fid,'%.3f pixels (photo %d over %d points)\n',mx,mxi,nPhoto(mxi));

fprintf(fid,[p2,'Point Precision\n']);

% Variance for each OP.
v=s.post.std.OP.^2;
% Mask fixed OPs.
v(~s.bundle.est.OP)=nan;
% Total variance.
tVar=sum(v,1);
% Total standard deviation.
tStd=sqrt(tVar);

fprintf(fid,[p3,'Total standard deviation (RMS of X/Y/Z std):\n']);

[mn,mni]=min(tStd);
[mx,mxi]=max(tStd);
fprintf(fid,[p4,'Minimum: ']);
fprintf(fid,'%.2g (OP %d)\n',mn,s.OP.id(mni));
fprintf(fid,[p4,'Maximum: ']);
fprintf(fid,'%.2g (OP %d)\n',mx,s.OP.id(mxi));

% Component-wise.
[mx,mxi]=max(sqrt(v),[],2);
for i=1:3
    fprintf(fid,[p3,'Maximum %c standard deviation: %.2g (OP %d)\n'],...
                 abs('X')+i-1,mx(i),s.OP.id(mxi(i)));
end

% Points with highest correlations.
fprintf(fid,[p3,'Points with high correlations\n']);
opCorr95=abs(vop)>0.95;
opCorr99=abs(vop)>0.99;
fprintf(fid,[p4,'Points with correlation above 95%%: %d\n'],nnz(opCorr95));
fprintf(fid,[p4,'Points with correlation above 99%%: %d\n'],nnz(opCorr99));
if nnz(opCorr95)
    fprintf(fid,[p4,'Points with highest correlations:\n']);
    [~,i]=sort(abs(vop),'descend');
    printed=[];
    j=1;
    while length(printed)<5 && j<=length(i)
        % Same OP may appear multiple times due to e.g. high X-Z,
        % Y-Z corr. Only print each OP once.
        if ~ismember(kop(i(j)),printed)
            printed(end+1)=kop(i(j));
            fprintf(fid,[p5,'Points %d: %.2f\n'], kop(i(j)),...
                    100*vop(i(j)));
        end
        j=j+1;
    end
end

fprintf(fid,[p2,'Point Angles\n']);

% Maximum angle between rays.
a=angles(s,'Computing angles')*180/pi;

fprintf(fid,[p3,'CP\n']);

% Ctrl pts with >0 rays
cpIx=find(s.prior.OP.isCtrl);
cpRays=sum(s.IP.vis(cpIx,:),2);
if any(cpRays==0)
    fprintf(fid,[p4,'Ignoring %d CP with 0 rays.\n'],nnz(cpRays==0));
end
% If we have any ctrl pts with rays.
if any(cpIx(cpRays>0))
    aCP=a(s.prior.OP.isCtrl);
    idCP=s.OP.id(s.prior.OP.isCtrl);
    labelCP=s.OP.label(s.prior.OP.isCtrl);
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
    fprintf(fid,[p4,'Minimum: %.1f degrees (CP %d%s)\n'],mn,idCP(mni),lmn);
    fprintf(fid,[p4,'Maximum: %.1f degrees (CP %d%s)\n'],mx,idCP(mxi),lmx);
    fprintf(fid,[p4,'Average: %.1f degrees\n'],NaNMean(aCP));
else
    fprintf(fid,[p4,'Minimum: -\n']);
    fprintf(fid,[p4,'Maximum: -\n']);
    fprintf(fid,[p4,'Average: -\n']);
end

fprintf(fid,[p3,'CCP\n']);
if any(s.prior.OP.isCheck)
    aCCP=a(s.prior.OP.isCheck);
    idCCP=s.OP.id(s.prior.OP.isCheck);
    labelCCP=s.OP.label(s.prior.OP.isCheck);
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
    fprintf(fid,[p4,'Minimum: %.1f degrees (CCP %d%s)\n'],mn,idCCP(mni),lmn);
    fprintf(fid,[p4,'Maximum: %.1f degrees (CCP %d%s)\n'],mx,idCCP(mxi),lmx);
    fprintf(fid,[p4,'Average: %.1f degrees\n'],mean(aCCP));
else
    fprintf(fid,[p4,'Minimum: -\n']);
    fprintf(fid,[p4,'Maximum: -\n']);
    fprintf(fid,[p4,'Average: -\n']);
end

fprintf(fid,[p3,'OP\n']);
isOP=~s.prior.OP.isCtrl & ~s.prior.OP.isCheck;
if any(isOP)
    aOP=a(isOP);
    idOP=s.OP.id(isOP);
    [mn,mni]=min(aOP);
    [mx,mxi]=max(aOP);
    fprintf(fid,[p4,'Minimum: %.1f degrees (OP %d)\n'],mn,idOP(mni));
    fprintf(fid,[p4,'Maximum: %.1f degrees (OP %d)\n'],mx,idOP(mxi));
    fprintf(fid,[p4,'Average: %.1f degrees\n'],mean(aOP));
    fprintf(fid,[p4,'Smallest angles (ID, angle [deg], vis in ' ...
                 'cameras)\n']);
    [ang,i]=sort(aOP);
    % Get all points with at least third min angle + 10% + 0.1 deg...
    angLimit=ang(min(3,end))*1.1+0.1;
    % ...but no point above 80 degress.
    angLimit=min(angLimit,80);
    % At least 3 points but obviously not more than all.
    nPts=min(max(nnz(ang<angLimit),3),length(ang));
    for j=1:nPts
        camVis=find(s.IP.vis(i(j),:));
        str=sprintf('%4d ',camVis);
        fprintf(fid,[p5,'%6d: %5.2f (%s)\n'],idOP(i(j)),ang(j),str(1:end-1));
    end
else
    fprintf(fid,[p4,'Minimum: -\n']);
    fprintf(fid,[p4,'Maximum: -\n']);
    fprintf(fid,[p4,'Average: -\n']);
end

fprintf(fid,[p2,'Ctrl measurements\n']);
if ~any(s.prior.OP.isCtrl)
    fprintf(fid,[p3,'none\n']);
else
    cIx=find(s.prior.OP.isCtrl);
    CPid=s.OP.id(cIx);

    fprintf(fid,[p3,'Prior\n']);
    fprintf(fid,[p3,'%6s, %8s, %8s, %8s, %8s, %8s, %8s, %s\n'],...
            'id','x','y','z','stdx','stdy','stdz','label');
    pos0=s.prior.OP.val(:,cIx);
    std0=s.prior.OP.std(:,cIx);
    for i=1:length(cIx)
        fprintf(fid,[p3,'%6d, %8.3f, %8.3f, %8.3f, %8.3g, ' ...
                     '%8.3g, %8.3g, %s\n'],CPid(i),pos0(:,i),std0(:,i),...
                s.OP.label{cIx(i)});
    end
    fprintf(fid,[p3,'Posterior\n']);
    fprintf(fid,[p3,'%6s, %8s, %8s, %8s, %8s, %8s, %8s, %4s, %s\n'],...
            'id','x','y','z','stdx','stdy','stdz','rays','label');
    pos1=s.OP.val(:,cIx);
    std1=s.post.std.OP(:,cIx);
    for i=1:length(cIx)
        fprintf(fid,[p3,'%6d, %8.3f, %8.3f, %8.3f, %8.3g, ' ...
                     '%8.3g, %8.3g, %4d, %s\n'],CPid(i),pos1(:,i),std1(:,i),...
                full(sum(s.IP.vis(cIx(i),:))),s.OP.label{cIx(i)});
    end
    fprintf(fid,[p3,'Diff (pos=abs diff, std=rel diff)\n']);
    fprintf(fid,[p3,'%6s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %4s, %s\n'],...
            'id','x','y','z','xy','xyz','stdx','stdy','stdz','rays','label');
    posd=pos1-pos0;
    stdd=((std1+eps)./(std0+eps)-1)*100;
    for i=1:length(cIx)
        fprintf(fid,[p3,'%6d, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %7.1f%%, ' ...
        '%7.1f%%, %7.1f%%, %4d, %s\n'],CPid(i),posd(:,i),norm(posd(1:2,i)),...
                norm(posd(:,i)),stdd(:,i),...
                full(sum(s.IP.vis(cIx(i),:))),s.OP.label{cIx(i)});
    end

    fprintf(fid,[p3,'Ctrl point delta\n']);
    diffNorm=sqrt(sum(posd.^2));
    [mx,i]=max(diffNorm);
    maxId=CPid(i);
    maxLabel=s.OP.label{cIx(i)};
    if ~isempty(maxLabel)
        maxLabel=[maxLabel,', '];
    end
    fprintf(fid,[p4,'Max: %.3f ou (%spt %d)\n'],mx,maxLabel,maxId);
    fprintf(fid,[p4,'Max X,Y,Z\n']);
    for i=1:3
        [mx,j]=max(abs(posd(i,:)));
        ll=s.OP.label{cIx(j)};
        if ~isempty(ll)
            ll=[ll,', '];
        end
        fprintf(fid,[p5,'%c: %.3f ou (%spt %d)\n'],'X'-1+i,mx,ll,CPid(j));
    end
    fprintf(fid,[p4,'RMS: %.3f ou (from %d items)\n'],...
            sqrt(mean(diffNorm.^2)),length(diffNorm));
end

fprintf(fid,[p2,'Check measurements\n']);
if ~any(s.prior.OP.isCheck)
    fprintf(fid,[p3,'none\n']);
else
    cIx=find(s.prior.OP.isCheck);
    CPid=s.OP.id(cIx);

    fprintf(fid,[p3,'Prior\n']);
    fprintf(fid,[p3,'%6s, %8s, %8s, %8s, %8s, %8s, %8s, %s\n'],...
            'id','x','y','z','stdx','stdy','stdz','label');
    pos0=s.prior.OP.val(:,cIx);
    std0=s.prior.OP.std(:,cIx);
    for i=1:length(cIx)
        fprintf(fid,[p3,'%6d, %8.3f, %8.3f, %8.3f, %8.3g, ' ...
                     '%8.3g, %8.3g, %s\n'],CPid(i),pos0(:,i),std0(:,i),...
                s.OP.label{cIx(i)});
    end
    fprintf(fid,[p3,'Posterior\n']);
    fprintf(fid,[p3,'%6s, %8s, %8s, %8s, %8s, %8s, %8s, %4s, %s\n'],...
            'id','x','y','z','stdx','stdy','stdz','rays','label');
    pos1=s.OP.val(:,cIx);
    std1=s.post.std.OP(:,cIx);
    for i=1:length(cIx)
        fprintf(fid,[p3,'%6d, %8.3f, %8.3f, %8.3f, %8.3g, ' ...
                     '%8.3g, %8.3g, %4d, %s\n'],CPid(i),pos1(:,i),std1(:,i),...
                full(sum(s.IP.vis(cIx(i),:))),s.OP.label{cIx(i)});
    end
    fprintf(fid,[p3,'Diff (pos=abs diff, std=rel diff)\n']);
    fprintf(fid,[p3,'%6s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %4s, %s\n'],...
            'id','x','y','z','xy','xyz','stdx','stdy','stdz','rays','label');
    posd=pos1-pos0;
    stdd=((std1+eps)./(std0+eps)-1)*100;
    for i=1:length(cIx)
        fprintf(fid,[p3,'%6d, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %7.1f%%, ' ...
        '%7.1f%%, %7.1f%%, %4d, %s\n'],CPid(i),posd(:,i),norm(posd(1:2,i)),...
                norm(posd(:,i)),stdd(:,i),...
                full(sum(s.IP.vis(cIx(i),:))),s.OP.label{cIx(i)});
    end

    fprintf(fid,[p3,'Check point delta\n']);
    diffNorm=sqrt(sum(posd.^2));
    [mx,i]=max(diffNorm);
    maxId=CPid(i);
    maxLabel=s.OP.label{cIx(i)};
    fprintf(fid,[p4,'Max: %.3f ou (%s, pt %d)\n'],mx,maxLabel,maxId);
    fprintf(fid,[p4,'Max X,Y,Z\n']);
    for i=1:3
        [mx,j]=max(abs(posd(i,:)));
        fprintf(fid,[p5,'%c: %.3f ou (%s, pt %d)\n'],'X'-1+i,mx,...
                s.OP.label{cIx(j)},CPid(j));
    end
    fprintf(fid,[p4,'RMS: %.3f ou (from %d items)\n'],...
            sqrt(mean(diffNorm.^2)),length(diffNorm));
end

fprintf(fid,'End of result file\n');

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


function m=NaNMean(v)
% Compute mean of non-nan values of v.

i=isnan(v);
if nnz(~i)==0
    m=nan;
else
    m=mean(v(~i));
end
