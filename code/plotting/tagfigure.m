function hh=tagfigure(tag,name)
%TAGFIGURE Create figure window based on tag name.
%
%   TAGFIGURE('tag') creates a figure with tag 'tag' unless such a figure
%   already exists.  In either case, the figure with tag 'tag' is brought
%   forward.  This correspond to the behaviuor of the FIGURE command,
%   except a string is used as the identifier.
%
%   TAGFIGURE('tag',TRUE) will also put the tag in the figure name.
%
%   H=TAGFIGURE(...) returns a handle to the figure.


narginchk(1,2);

if nargin<2, name=false; end

% Search for an already opened figure window.
h=findobj('type','figure','tag',tag);
if isempty(h)
    % Create one if it didn't exist.
    h=figure('tag',tag);
else
    figure(h)
end

if name
    set(h,'name',tag);
end

if nargout>0
    hh=h;
end
