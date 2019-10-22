function [ok,msg]=checkversionattr(attr,fieldName,lastKnownVer,firstKnownVer)
%CHECKVERSIONATTR Check version attribute is present and known
%
%   OK=CHECKVERSIONATTR(ATTR,FIELDNAME,LASTVER) checks that the field
%   FIELDNAME is present in the structure ATTR and that the version
%   string therein is compared to the version string LASTVER. The
%   function returns OK as TRUE if the version is less than or equal
%   to LASTVER. Use CHECKVERSIONATTR(ATTR,FIELDNAME,LASTVER,FIRSTVER)
%   to also check that the version is more than or equal to
%   FIRSTVER. Use a blank version string for LASTVER or FIRSTVER to
%   defer the corresponding test.
%
%   Note: The version strings can contain at most four components,
%   e.g., '1.6.8.1645'.

if nargin<4, firstKnownVer=''; end

[ok,msg]=checkxmlfields(attr,fieldName);

if ~ok
    return
end

% Parse current version strings.
curVer=attr.(fieldName);
curParts=sscanf(curVer,'%d.%d.%d.%d');
if length(curParts)<4, curParts(4)=0; end

% Parse limiting strings if specified.
if ~isempty(lastKnownVer)
    lastKnownParts=sscanf(lastKnownVer,'%d.%d.%d.%d');
    if length(lastKnownParts)<4, lastKnownParts(4)=0; end
else
    lastKnownParts=+inf;
end

if ~isempty(firstKnownVer)
    firstKnownParts=sscanf(firstKnownVer,'%d.%d.%d.%d');
    if length(firstKnownParts)<4, firstKnownParts(4)=0; end
else
    firstKnownParts=-inf;
end

msg='';
ok=true;

% Compare versions
for i=1:4
    if curParts(i)>lastKnownParts(i)
        ok=false;
        msg=sprintf('Unknown version %s (last known version is %s)',...
                    curVer,lastKnownVer);
        return;
    elseif curParts(i)<lastKnownParts(i)
        return;
    else
        % Version parts are equal. Continue and check the next part.
    end
end

for i=1:4
    if curParts(i)<firstKnownParts(i)
        ok=false;
        msg=sprintf('Unknown version %s (first known version is %s)',...
                    curVer,firstKnownVer);
        return;
    elseif curParts(i)>lastKnownParts(i)
        return;
    else
        % Version parts are equal. Continue and check the next part.
    end
end
