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

p=[];

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
            groups=doc.getChildNodes;
            for iGroup=0:groups.getLength()-1
                group=groups.item(iGroup);
                switch char(group.getNodeName())
                  case 'group'
                    disp('Group');
                  case 'report'
                    disp('Report');
                  case '#text'
                    % Ignore #text as it is usually just linebreaks.
                  otherwise
                    warning(sprintf('Expected group or report, got %s.',...
                                  char(group.getNodeName())));
                end
            end
            
        end
    end
end

% % Loop over all elements.
% for d=0:doc.getLength-1
%     i=doc.item(d);

 
% %Note that the item list index is zero-based.
% for i=0:allListItems.getLength-1
%     thisListItem=allListItems.item(i);
%     childNode = thisListItem.getFirstChild;
 
%     while ~isempty(childNode)
%         %Filter out text, comments, and processing instructions.
%         if childNode.getNodeType == childNode.ELEMENT_NODE
%             %Assume that each element has a single org.w3c.dom.Text child
%             childText = char(childNode.getFirstChild.getData);
%             switch char(childNode.getTagName)
%               case 'label' ; itemFound = strcmp(childText,infoLabel);
%               case 'callback' ; infoCbk = childText;
%             end
%         end
%         childNode = childNode.getNextSibling;
%     end
%     if itemFound break; else infoCbk = ''; end
% end
% disp(sprintf('Item "%s" has a callback of "%s".',infoLabel,infoCbk))
 