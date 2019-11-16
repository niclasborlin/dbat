function s=parseoutput(s,output)
%PARSEOUTPUT Parse and execute output section of a DBAT XML file.
%
%    S=PARSEOUTPUT(S,OUTPUT) parses the OUTPUT XML block from a DBAT
%    XML script file and generates the requested plots and files.
%
%See also: PARSEINPUT, PARSEOPS, RUNDBATSCRIPT.
    
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

if isfield(output,'files') && isfield(output.files,'file')
    files=output.files;
    if ~iscell(files)
        % Single file
        files={files};
    end
else
    files={};
end

% For each plot...
for i=1:length(plots)
    plot=plots{i};
    if ~isfield(plot,'Text')
        error('DBAT XML script operations error: Malformed plot no %d',i);
    end
    switch plot.Text
      case 'image'
        PlotImage(s,plot)
      case 'image_stats'
        warning('%s plot not implemented yet',plot.Text)
      case 'op_stats'
        warning('%s plot not implemented yet',plot.Text)
      case 'coverage'
        warning('%s plot not implemented yet',plot.Text)
      case 'params'
        warning('%s plot not implemented yet',plot.Text)
      case 'iteration_trace'
        warning('%s plot not implemented yet',plot.Text)
      otherwise
        error('DBAT XML script operations error: Unknown operation %s',...
              plot.Text)
    end
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
