function [rr,s0,psz]=stpierrebundledemo_ps
%STPIERREBUNDLEDEMO_PS StPierre PhotoScan bundle demo for DBAT.
%
%   STPIERREBUNDLEDEMO_PS runs the bundle on the PhotoScan project
%   file for the STPIERRE data set, used in [1]. The PhotoScan EO and
%   OP values are used as initial values. The datum is defined via
%   weighted control points. The Gauss-Newton-Armijo damping scheme
%   of [2] is used for the bundle.
%
%   References:
%       [1] A. Murtiyoso and P. Grussenmeyer and N. Börlin (2017).
%           "Reprocessing Close Range Terrestrial and UAV
%           Photogrammetric Projects with the DBAT Toolbox for
%           Independent Verification and Quality Control",
%           International Archives of the Photogrammetry, Remote
%           Sensing and Spatial Information Sciences, XLII-2/W8, p.
%           171-177. Paper presented at the LowCost 3D Workshop,
%           Hamburg 28-29 November 2017.
%       [2] Börlin and Grussenmeyer (2013). "Bundle adjustment with
%           and without damping", Photogrammetric Record, vol.
%           28(144):396-415.

damping='gna';

% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

% Base dir with input files.
inputDir=fullfile(curDir,'data','hamburg2017','stpierre');

% Name of control point file.
cpName=fullfile(inputDir,'ctrl_StPierre_weighted.txt');

% PhotoModeler text export file and report file.
inputFile=fullfile(inputDir,'psprojects','C5.psz');
% Report file name.
reportFile=fullfile(inputDir,'dbatexports','stpierrePS_C5-dbatreport.txt');

fprintf('Loading PhotoScan project file %s...',inputFile);
psz=loadpsz(inputFile);
fprintf('done.\n');

[prob,pmReport,pts3d,pts2d]=ps2pmstruct(psz);

probRaw=prob; %#ok<NASGU>
if any(isnan(cat(2,prob.images.imSz)))
    error('Image sizes unknown!');
end
disp('done.')

markPtsStd=unique(prob.markPts(:,5:6));
if isscalar(markPtsStd)
    disp('Same stdev on all mark points:')
    disp(markPtsStd);
else
    disp('Different stdev on some mark points:')
    disp(markPtsStd);
end    
%fprintf('Setting all mark pts std to 1.\n');
%prob.markPts(:,5:6)=1;

fprintf('Loading control point file %s...',cpName);
ctrlPts=loadcpt(cpName);
fprintf('done.\n');

% Convert loaded PhotoModeler data to DBAT struct.
s0=prob2dbatstruct(prob);
ss0=s0;

[i1,j1,i2,j2,j3]=matchcpt(s0,ctrlPts,'label');

% Compute rigid-body transformation for control points.
TCP=rigidalign(s0.prior.OP.val(:,i1),ctrlPts.pos(:,j1));
if length(i2)>=3
    TCC=rigidalign(s0.prior.OP.val(:,i2),ctrlPts.pos(:,j2));
else
    TCC=[];
end

% Display info about CPT
fprintf('\nFound %d control points:\n',length(i1))
disp(ctrlPts.pos(:,j1))
disp(ctrlPts.name(j1))
angCP=acosd(min(max((trace(TCP(1:3,1:3))-1)/2,0),1));
fprintf('Rotation %.1f degress, translation %.1f object units.\n',...
        angCP,norm(TCP(1:3,4)));
fprintf('Max abs pos diff=%g\n',max(max(abs(s0.prior.OP.val(:,i1)-ctrlPts.pos(:,j1)))));
fprintf('Max rel std diff=%.1f%%\n',(exp(max(max(abs(log(s0.prior.OP.val(:,i1))-log(ctrlPts.pos(:,j1))))))-1)*100);

% Display info about CC
fprintf('\nFound %d check points:\n',length(i2))
disp(ctrlPts.pos(:,j2))
disp(ctrlPts.name(j2))
if ~isempty(TCC)
    angCP=acosd(min(max((trace(TCC(1:3,1:3))-1)/2,0),1));
    fprintf('Rotation %.1f degress, translation %.1f object units.\n',...
            angCP,norm(TCC(1:3,4)));
end
fprintf('Max abs pos diff=%g\n',max(max(abs(s0.prior.OP.val(:,i2)-ctrlPts.pos(:,j2)))));
fprintf('Max rel std diff=%.1f%%\n',(exp(max(max(abs(log(s0.prior.OP.val(:,i2))-log(ctrlPts.pos(:,j2))))))-1)*100);

% Switch to Forward/Computer Vision lens distortion model for all cameras.
s0.IO.model.distModel(:)=-1;

% Self-calibration of all camera parameters except skew.
s0=setcamvals(s0,'loaded');
s0=setcamest(s0,'all','not','sk','as');

% Use supplied EO data as initial values.

% Use estimated OP values as initial.

% Datum is given by CPs.
% Do nothing.

fprintf('\nRunning the bundle with damping %s...\n',damping);

% Run the bundle.
[result,ok,iters,sigma0,E]=bundle(s0,damping,'trace');
    
if ok
    fprintf('Bundle ok after %d iterations with sigma0=%.2f (%.2f pixels)\n',...
            iters,sigma0,sigma0*s0.IP.sigmas(1));
else
    fprintf(['Bundle failed after %d iterations (code=%d). Last sigma0 estimate=%.2f ' ...
             '(%.2f pixels)\n'],iters,E.code,sigma0,sigma0*s0.IP.sigmas(1));
end

% Pre-factorize posterior covariance matrix for speed.
E=bundle_cov(result,E,'prepare');

[COP,result]=bundle_result_file(result,E,reportFile);

fprintf('\nBundle result file %s generated.\n',reportFile);

% Don't plot iteration history for the 26000+ object points.
h=plotparams(result,E,'noop');

h=plotcoverage(result,true);

h=plotimagestats(result,E);

h=plotopstats(result,E,COP);

fig=tagfigure(sprintf('damping=%s',damping));
fprintf('Displaying bundle iteration playback for method %s in figure %d.\n',...
        E.damping.name,double(fig));
plotnetwork(result,E,'trans','up','align',1,'title',...
            ['Damping: ',damping,'. Iteration %d of %d'], ...
            'axes',fig);

if nargout>0
    rr=result;
end

imName='';
imNo=1;
% Check if image files exist.
if exist(s0.proj.imDir,'dir')
    % Handle both original-case and lower-case file names.
    imNames={s0.EO.name{imNo},lower(s0.EO.name{imNo}),upper(s0.EO.name{imNo})};    
    imNames=fullfile(s0.proj.imDir,imNames);
    imExist=cellfun(@(x)exist(x,'file')==2,imNames);
    if any(imExist)
        imName=imNames{find(imExist,1,'first')};
    end
else
    warning('Image directory %s does not exist.',s0.proj.imDir);
end

if exist(imName,'file')
    fprintf('Plotting measurements on image %d.\n',imNo);
    imFig=tagfigure('image');
    h=[h;imshow(imName,'parent',gca(imFig))];
    pts=s0.IP.val(:,s0.IP.ix(s0.IP.vis(:,imNo),imNo));
    ptsId=s0.OP.id(s0.IP.vis(:,imNo));
    isCtrl=s0.prior.OP.isCtrl(s0.IP.vis(:,imNo));
    % Plot non-control points as red crosses.
    if any(~isCtrl)
        line(pts(1,~isCtrl),pts(2,~isCtrl),'marker','x','color','r',...
             'linestyle','none','parent',gca(imFig));
    end
    % Plot control points as black-yellow triangles.
    if any(isCtrl)
        line(pts(1,isCtrl),pts(2,isCtrl),'marker','^','color','k',...
             'markersize',2,'linestyle','none','parent',gca(imFig));
        line(pts(1,isCtrl),pts(2,isCtrl),'marker','^','color','y',...
             'markersize',6,'linestyle','none','parent',gca(imFig));
    end
    for i=1:length(ptsId)
        text(pts(1,i),pts(2,i),int2str(ptsId(i)),'horizontal','center',...
             'vertical','bottom','color','b','parent',gca(imFig));
    end
end
