function [data,err]=loadtextfile(filename,format,headerChar)
%LOADTEXTFILE Load data in text file with headers.
%
%   DATA=LOADTEXTFILE(FILENAME,FORMAT) uses TEXTSCAN to load data in
%   the text file FILENAME with format FORMAT. Header lines, i.e.,
%   blank lines or lines with a hash (#) character as the first
%   non-whitespace character are ignored. Data is returned from
%   TEXTSCAN without any processing.
%
%   DATA=LOADTEXTFILE(FILENAME,FORMAT,HEADERCHAR) uses HEADERCHAR
%   instead of #.
%
%   [DATA,ERR]=LOADTEXTFILE(...) will return an error string ERR on
%   error instead of throwing an error.
%
%See also: TEXTSCAN.

if nargin<3, headerChar='#'; end
headerChar=headerChar(1);

data=[];
err='';

% Open file
[fid,msg]=fopen(filename,'rt');
if fid<0
    err=sprintf('Failed to open file %s: %s',filename,msg);
    if nargout<2, error(err); else, return; end %#ok<SPERR>
end

% Read lines and check if they are headers.
fpos=0;

while ~feof(fid)
    % Remember position before reading next line.
    fpos=ftell(fid);
    if fpos<0
        err=sprintf('ftell failed on open file %s',filename);
        fclose(fid);
        if nargout<2, error(err); else, return; end %#ok<SPERR>
    end
    % Read next line.
    line=fgetl(fid);
    if isempty(line)
        % Check error condition.
        [msg,errNum]=ferror(fid);
        if isempty(msg)
            % No error, just a blank line. Read past it.
            continue;
        else
            err=sprintf('Error reading line from %s: %d %s',filename,errNum,msg);
            fclose(fid);
            if nargout<2, error(err); else, return; end %#ok<SPERR>
        end
    end
    i=find(~isspace(line),1);
    if ~isempty(i) && line(i)==headerChar
        % Header line, read next line.
        continue;
    else
        % Not a header line, we're done.
        break;
    end
end

% Reset position to before last line.
status=fseek(fid,fpos,'bof');
if status<0
    % Check error condition.
    [msg,errNum]=ferror(fid);
    if ~isempty(msg)
        err=sprintf('Error seeking in %s: %d %s',filename,errNum,msg);
        fclose(fid);
        if nargout<2, error(err); else, return; end %#ok<SPERR>
    end
end

data=textscan(fid,format);
fclose(fid);
