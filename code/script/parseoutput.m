function s=parseoutput(s,output,docFile)
%PARSEOUTPUT Parse and execute output section of a DBAT XML file.
%
%    S=PARSEOUTPUT(S,OUTPUT,DOCFILE) parses the OUTPUT XML block from
%    a DBAT XML script file and generates the requested plots and
%    files. The string DOCFILE should contain the path name of the XML
%    file and is used to determine base directories.
%
%See also: PARSEINPUT, PARSEOPS, PARSEOUTPUTFILES, RUNDBATSCRIPT.
    
narginchk(3,3);

% Known operations
knownFields={'plots','files','c'};

% All fields are optional.
[ok,msg]=checkxmlfields(output,knownFields,false(size(knownFields)));
if ~ok, error('DBAT XML script operations error: %s',msg); end

% Extract plots and files.
if isfield(output,'plots') && isfield(output.plots,'plot')
    plots=output.plots.plot;
    if ~iscell(plots)
        % Single plot
        plots={plots};
    end
else
    plots={};
end

% For each plot...
for i=1:length(plots)
    plot=plots{i};
    if ~isfield(plot,'Text')
        error('DBAT XML script operations error: Malformed plot no %d',i);
    end
    switch strip(plot.Text)
      case 'image'
        PlotImage(s,plot)
      case 'image_stats'
        PlotImageStats(s,plot)
      case 'op_stats'
        PlotOPStats(s,plot);
      case 'coverage'
        PlotCoverage(s,plot);
      case 'params'
        PlotParams(s,plot);
      case 'iteration_trace'
        PlotIterationTrace(s,plot);
      otherwise
        error('DBAT XML script operations error: Unknown operation %s',...
              plot.Text)
    end
end

% Parse and generate output files.
if isfield(output,'files')
    s=parseoutputfiles(s,output.files,docFile);
end


function PlotImage(s,op)
%Show image specified by the id attribute of op and plot
%observations on top.

if isfield(op,'Attributes') && isfield(op.Attributes,'id')
    imNo=sscanf(op.Attributes.id,'%d');
else
    imNo=nan;
end

if isnan(imNo)
    error('DBAT XML script plot/image error: Missing image id');
end

% Construct path name
imName=fullfile(s.proj.imDir,s.EO.name{imNo});
if exist(imName,'file')
    fprintf('Plotting measurements on image %d.\n',imNo);
    
    % Label as figure name and plot title. Full path as xlabel.
    imFig=tagfigure(sprintf('image-%d',imNo));
    imshow(imName,'parent',gca(imFig));
    h=gca(imFig);
    set(imFig,'name',sprintf('Image %d - %s',imNo,s.EO.label{imNo}));
    title(h,sprintf('Image %d - %s',imNo,s.EO.label{imNo}));
    xlabel(h,imName);
    
    % Extract points measured in this image.
    pts=s.IP.val(:,s.IP.ix(s.IP.vis(:,imNo),imNo));
    ptsId=s.OP.id(s.IP.vis(:,imNo));
    isCtrl=s.prior.OP.isCtrl(s.IP.vis(:,imNo));
    labels=s.OP.label(s.IP.vis(:,imNo));
    
    % Plot non-control points as red crosses.
    if any(~isCtrl)
        line(pts(1,~isCtrl),pts(2,~isCtrl),'marker','x','color','r',...
             'linestyle','none','parent',h);
    end

    % Plot control points as black-yellow triangles.
    if any(isCtrl)
        line(pts(1,isCtrl),pts(2,isCtrl),'marker','^','color','k',...
             'markersize',2,'linestyle','none','parent',gca(imFig));
        line(pts(1,isCtrl),pts(2,isCtrl),'marker','^','color','y',...
             'markersize',6,'linestyle','none','parent',gca(imFig));
    end
    
    % Plot point ids above points.
    for i=1:length(ptsId)
        text(pts(1,i),pts(2,i),int2str(ptsId(i)),'horizontal','center',...
             'vertical','bottom','color','b','parent',gca(imFig));
    end
    
    % Plot point labels below points.
    for i=find(cellfun(@length,labels))
        text(pts(1,i),pts(2,i),labels{i},'horizontal','center',...
             'vertical','top','color','c','parent',gca(imFig));
    end
end


function PlotParams(s,op)
%Show evolution of parameter estimates.

h=plotparams(s,s.bundle.info);


function PlotImageStats(s,op)
%Show image statistics.

if isempty(s.post.std.EO)
    CEO=bundle_cov(s,s.bundle.info,'CEO');
    s.post.std.EO=sqrt(reshape(diag(CEO),size(s.EO.val,1),[]));
    
    [mEO,nEO]=size(s.EO.val);
    s.post.cov.EO=nan(mEO,mEO,nEO);
    for i=1:nEO
        ix=(i-1)*mEO+(1:mEO);
        s.post.cov.EO(:,:,i)=CEO(ix,ix);
    end
end

h=plotimagestats(s,s.bundle.info);


function PlotOPStats(s,op)
% Plot object point statistics.

maxOP=inf;

if isfield(op,'Attributes')
    if isfield(op.Attributes,'max_op')
        maxOP=sscanf(op.Attributes.max_op,'%d');
    end
end

if size(s.OP.val,2)>maxOP
    warning('Number of OP (%d) exceed op_stats limit of %d. Not plotting.',...
            size(s.OP.val,2),maxOP);
else
    if isempty(s.post.std.OP)
        COP=bundle_cov(s,s.bundle.info,'COP');
        s.post.std.OP=sqrt(reshape(diag(COP),size(s.OP.val,1),[]));
    
        [mOP,nOP]=size(s.OP.val);
        s.post.cov.OP=nan(mOP,mOP,nOP);
        for i=1:nOP
            ix=(i-1)*mOP+(1:mOP);
            s.post.cov.OP(:,:,i)=COP(ix,ix);
        end
    end

    h=plotopstats(s,s.bundle.info);
end

function PlotIterationTrace(s,op)
% Plot iteration trace.

doPause='on';
camSize=1;

if isfield(op,'Attributes')
    if isfield(op.Attributes,'cam_size')
        camSize=sscanf(op.Attributes.cam_size,'%f');
    end
end

fig=tagfigure('trace');
damping=s.bundle.info.damping.name;
fprintf('Displaying bundle iteration playback for method %s in figure %d.\n',...
        damping,double(fig));
h=plotnetwork(s,s.bundle.info,'title',...
              ['Damping: ',damping,'. Iteration %d of %d'], ...
              'axes',fig,'pause',doPause,'camsize',camSize); 


function PlotCoverage(s,op)
% Plot image coverage.

convexHull=false;

if isfield(op,'Attributes')
    if isfield(op.Attributes,'convex_hull')
        convexHull=strcmp(strip(op.Attributes.convex_hull),'true');
    end
end

h=plotcoverage(s,convexHull);
