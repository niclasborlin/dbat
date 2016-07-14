function DOPRINTDEMOFIGURES(figs,names,dir)
%DOPRINTDEMOFIGURES Print demo figures to disk.
%
%   DOPRINTDEMOFIGURES(FIGS,NAMES), where FIGS is an N-vector with figure
%   and/or axes handles and NAMES is an N-cell array with .EPS file names,
%   writes the figures to files in the default manual illustration
%   directory. If FIGS contains only one handle, NAMES may be a single
%   string.
%
%   DOPRINTDEMOFIGURES(FIGS,NAMES,DIR) writes the files to the directory DIR
%   instead.

% $Id$

if nargin<3
    % Default directory.
    dir=fullfile(fileparts(mfilename('fullpath')),'..','..','doc',...
                 'manual','ill');
end

if isscalar(figs) && ~iscell(names) && ischar(names)
    names={names};
end

if length(figs)~=length(names)
    error('FIGS and NAMES must be of equal length');
end

for i=1:length(figs)
    % Get parent figure handle, if necessary. Add safety stop for root
    % handle (should never happen).
    h=figs(i);
    while ~strcmp(get(h,'type'),'figure') || h==0
        h=get(h,'parent');
    end
    if h~=0
        print(h(i),'-depsc2',fullfile(dir,names{i}));
    end
end
