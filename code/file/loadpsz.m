function [p,doc]=loadpsz(fName)
%LOADPSZ Load Photoscan 1.1.6 PSZ camera calibration file.
%
%

% Unpack zip archive into temporary folder.
tmpDir=tempname;
files=unzip(fName,tmpDir);

% Find doc.xml file.
docFile=files{strcmp(fullfile(tmpDir,'doc.xml'),files)};

% Read the doc.xml file.
obj=xmlread(docFile);

p=struct('chunks',[]);

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
            warning('No chunks found');
        else
            % Get document version.
            docAttrs=doc.getAttributes();
            docVersionStr=docAttrs.getNamedItem('version');
            docVersion=regexp(char(docVersionStr.toString()),'version="(.*)"','match','once');
            % Document should have chunks.
            chunks=doc.getChildNodes;
            for i=0:chunks.getLength()-1
                child=chunks.item(i);
                switch char(child.getNodeName())
                  case '#text'
                    % Ignore #text as it is usually just linebreaks.
                  case 'chunks'
                    c=ParseChunks(child);
                  otherwise
                    warning(sprintf('Expected chunks, got %s.',...
                                  char(child.getNodeName())));
                end
            end
        end
    end
end

function c=ParseChunks(chunks)
% Parse chunks.

c=struct('chunk',[]);

if ~chunks.hasChildNodes()
    warning('Chunks without children')
    return;
end

% Chunk should have photos as children.
children=chunks.getChildNodes;

for iChild=0:chunks.getLength()-1
    child=children.item(iChild);
    switch char(child.getNodeName())
      case '#text'
        % Ignore #text as it is usually just linebreaks.
      case 'chunk'
        cc=ParseChunk(child);
        if isempty(c.chunk)
            c.chunk=cc;
        else
            c.chunk(end+1)=cc;
        end
      otherwise
        warning(sprintf('Unexpected chunk child %s.',...
                        char(child.getNodeName())));
    end
end

function c=ParseChunk(chunk)
% Parse chunks.

c=struct('chunk',[]);

if ~chunk.hasChildNodes()
    warning('Chunk without children')
    return;
end

% Chunk should have photos as children.
children=chunk.getChildNodes;

for iChild=0:chunk.getLength()-1
    child=children.item(iChild);
    name=char(child.getNodeName());
    switch name
      case '#text'
        % Ignore #text as it is usually just linebreaks.
      case 'sensors'
        cc=ParseSensors(child);
      case 'cameras'
      case 'markers'
      case 'frames'
      case 'transform'
      case 'region'
      case 'settings'
      case 'meta'
        cc=ParseChunk(child);
        if isempty(c.chunk)
            c.chunk=cc;
        else
            c.chunk(end+1)=cc;
        end
      otherwise
        warning(sprintf('Unexpected chunk child %s.',...
                        char(child.getNodeName())));
    end
end
