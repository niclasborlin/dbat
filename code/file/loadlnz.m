function [p,doc]=loadlnz(fName)
%LOADLNZ Load Photoscan LNZ camera calibration file.
%
%

% Unpack zip archive into temporary folder.
tmpDir=tempname;
files=unzip(fName,tmpDir);
if length(files)~=1
    error('Expected single .xml file.');
end

% Read the doc.xml file.
obj=xmlread(files{1});

% Clean up the file and dir.
delete(files{1});
rmdir(tmpDir);

p=struct('group',[],'report',[]);

if ~obj.hasChildNodes
    warning('No document found');
else
    % File should have at least one document.
    docs=obj.getChildNodes;
    for iDoc=0:docs.getLength()-1
        doc=docs.item(iDoc);
        if ~strcmp(doc.getNodeName(),'document')
            warning('Not a document');
        elseif ~obj.hasChildNodes
            warning('No groups/reports found');
        else
            % Document should have groups and reports.
            docChildren=doc.getChildNodes;
            for iChild=0:docChildren.getLength()-1
                child=docChildren.item(iChild);
                switch char(child.getNodeName())
                  case 'group'
                    group=ParseGroup(child);
                    if isempty(p.group)
                        p.group=group;
                    else
                        p.group(end+1)=group;
                    end
                  case 'report'
                    disp('Report');
                    report=ParseReport(child);
                    if isempty(p.report)
                        p.report=report;
                    else
                        p.report(end+1)=report;
                    end
                  case '#text'
                    % Ignore #text as it is usually just linebreaks.
                  otherwise
                    warning(sprintf('Expected group or report, got %s.',...
                                  char(child.getNodeName())));
                end
            end
        end
    end
end

% --- Internal function ParseGroup ---
function g=ParseGroup(group)
% Parse group.

g=struct('photo',[]);

if ~group.hasChildNodes()
    warning('Group without children')
    return;
end

% Group should have photos as children.
children=group.getChildNodes;

for iChild=0:group.getLength()-1
    child=children.item(iChild);
    switch char(child.getNodeName())
      case 'photo'
        photo=ParsePhoto(child);
        if isempty(g.photo)
            g.photo=photo;
        else
            g.photo(end+1)=photo;
        end
      case '#text'
        % Ignore #text as it is usually just linebreaks.
      otherwise
        warning(sprintf('Unexpected group child %s.',...
                        char(child.getNodeName())));
    end
end
            

% --- Internal function ParsePhoto ---
function p=ParsePhoto(photo)

p=[];

if ~photo.hasChildNodes()
    warning('Photo without children')
    return;
end

% Photos should have lots of children.
children=photo.getChildNodes;

% Create blank structure with 1000 preallocated points.
p=struct('location','','info',struct,'corner',zeros(1000,5),'transform',[],'calibration',[]);

nCorner=0;
for iChild=0:photo.getLength()-1
    child=children.item(iChild);
    switch char(child.getNodeName())
      case 'location'
        a=child.getAttributes;
        % Extract image path.
        p.location=char(a.getNamedItem('path').getValue);
      case 'meta'
        a=child.getAttributes;
        % Extract name and value for metadata.
        p.info=setfield(p.info,char(a.getNamedItem('name').getValue),...
                               char(a.getNamedItem('value').getValue));
      case 'corner'
        % Extract image and object space coordinates and status for the corner.
        nCorner=nCorner+1;
        a=child.getAttributes();
        imgxy=[sscanf(char(a.getNamedItem('img_x').getValue),'%g'),...
               sscanf(char(a.getNamedItem('img_y').getValue),'%g')];
        objxy=[sscanf(char(a.getNamedItem('obj_x').getValue),'%g'),...
               sscanf(char(a.getNamedItem('obj_y').getValue),'%g')];
        valid=strcmp(char(a.getNamedItem('valid').getValue),'true');
        % If preallocated space is full, make a block expansion.
        if size(p.corner,1)<nCorner
            p.corner(end+1000,1)=0;
        end
        p.corner(nCorner,:)=[imgxy,objxy,valid];
      case 'transform'
        p.transform=ParseTransform(child);
      case 'calibration'
        p.calibration=ParseCalibration(child);
      case '#text'
        % Ignore #text as it is usually just linebreaks.
      otherwise
        warning(sprintf('Unexpected photo child %s.',...
                        char(child.getNodeName())));
    end
end

p.corner=p.corner(1:nCorner,:);


% --- Internal function ParseTransform ---
function t=ParseTransform(transform)

% Transformation is 4-by-4 matrix, row-wise.
t=reshape(sscanf(char(transform.getFirstChild.getTextContent()),'%g'),4,4)';


% -- Internal function ParseCalibration ---
function c=ParseCalibration(calibration)

if ~calibration.hasChildNodes()
    warning('Calibration without children')
    return;
end

% Calibrations should have lots of children.
children=calibration.getChildNodes;

c=struct('width',nan,'height',nan,'fx',nan,'fy',nan,'cx',nan,'cy',nan,...
         'skew',0,'k1',0,'k2',0,'k3',0,'k4',0,'p1',0,'p2',0);

for iChild=0:calibration.getLength()-1
    child=children.item(iChild);
    switch char(child.getNodeName())
      case {'width','height','fx','fy','cx','cy','skew','k1','k2','k3','k4','p1','p2'}
        val=sscanf(char(child.getFirstChild.getTextContent()),'%g');
        c=setfield(c,char(child.getNodeName()),val);
      case '#text'
        % Ignore #text as it is usually just linebreaks.
      otherwise
        warning(sprintf('Unexpected calibration child node %s.',...
                        char(child.getNodeName())));
    end
end


function r=ParseReport(report)

r=[];

if ~report.hasChildNodes()
    warning('Report without children')
    return;
end

% Report should have lots of children.
children=report.getChildNodes;

c=struct('width',nan,'height',nan,'focal',nan,'fx',nan,'fy',nan,...
         'cx',nan,'cy',nan,'skew',0,'k1',0,'k2',0,'k3',0,'k4',0,'p1',0,'p2',0);

for iChild=0:children.getLength()-1
    child=children.item(iChild);
    switch char(child.getNodeName())
      case 'record'
        record=child;
        a=record.getAttributes();
        c.width=sscanf(char(a.getNamedItem('width').getValue),'%g');
        c.height=sscanf(char(a.getNamedItem('height').getValue),'%g');
        c.focal=sscanf(char(a.getNamedItem('flength').getValue),'%g');
        recChildren=record.getChildNodes();
        for iRecChild=0:recChildren.getLength()-1
            recChild=recChildren.item(iRecChild);
            a=recChild.getAttributes();
            switch char(recChild.getNodeName())
              case {'fx','fy','cx','cy','skew','k1','k2','k3','k4','p1','p2'}
                val=sscanf(char(a.getNamedItem('value').getValue),'%g');
                err=sscanf(char(a.getNamedItem('error').getValue),'%g');
                c=setfield(c,char(recChild.getNodeName()),...
                             struct('value',val,'error',err));
              case '#text'
                % Ignore #text as it is usually just linebreaks.
              otherwise
                warning(sprintf('Unexpected record child node %s.',...
                                char(child.getNodeName())));
            end
        end
      case '#text'
        % Ignore #text as it is usually just linebreaks.
      otherwise
        warning(sprintf('Unexpected calibration child node %s.',...
                        char(child.getNodeName())));
    end
end

r.record=c;