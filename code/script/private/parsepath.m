function p=parsepath(p,baseDir)
%PARSEPATH Parse a DBAT file path and return an absolute path
%
%   P=PARSEPATH(P,BASEDIR) prepends BASEDIR to P is P is a relative
%   path, otherwise returns P unchanged. P is an absolute path if it
%   starts with slash, backslash, or 'X:', where X is any letter.
%
%   Any leading and trailing whitespace is stripped from P before
%   processing.

narginchk(2,2)

% Remove all leading and trailing whitespace.
p=strip(p);

if isempty(p) || p(1)=='/' || p(1)=='\' || ...
        (length(p)>=2 && isstrprop(p(1),'alpha') && p(2)==':')
    % Path is absolute, do nothing.
    return
end

p=fullfile(baseDir,p);

