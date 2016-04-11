function comparepsz(psz1,psz2,verb)
%COMPAREPSZ Compare two .psz files.
%
%    COMPAREPSZ(PSZ1,PSZ2) compare the content of the .psz files PSZ1
%    and PSZ2, printing the changes. The content of the unpacked
%    directories will be overwritten.

if nargin<3, verb=false; end

unpackpsz(psz1)
unpackpsz(psz2)

pszDir1=fileparts(psz1);
pszDir2=fileparts(psz2);

pszUnpackedDir1=fullfile(pszDir1,'unpacked');
pszUnpackedDir2=fullfile(pszDir2,'unpacked');

docFile1=fullfile(pszUnpackedDir1,'doc.xml');
docFile2=fullfile(pszUnpackedDir2,'doc.xml');

s1=xml2struct2(docFile1);
s2=xml2struct2(docFile2);

eq=compare(s1,s2,'s',verb);


function eq=compare(s1,s2,level,verb)

if isequalwithequalnans(s1,s2)
    if verb
        fprintf('%s: equal.\n',level);
    end
    eq=true;
    return
end

numericStrings={'Text','value','x','y','z','sx','sy','sz','sxy','sxyz'};

if isstruct(s1) && isstruct(s2)
    fld1=fields(s1);
    fld2=fields(s2);
    flds=unique({fld1{:},fld2{:}});
    for i=1:length(flds)
        if xor(isfield(s1,flds{i}),isfield(s2,flds{i}))
            eq=false;
            if isfield(s1,flds{i})
                fprintf('Field %s.%s missing in s2.\n',level,flds{i});
            else
                fprintf('Field %s.%s missing in s1.\n',level,flds{i});
            end
        elseif strcmp(flds{i},'Text')
            % Special treatment of .Text fields.
            v1=sscanf(getfield(s1,flds{i}),'%g ');
            v2=sscanf(getfield(s2,flds{i}),'%g ');
            if length(v1)~=length(v2)
                fprintf('%s.Text: length difference (%d vs. %d).\n',level,...
                        length(v1),length(v2));
                eq=false;
            elseif isempty(getfield(s1,flds{i}))
                % Silent return.
                eq=true;
            else
                df=max(abs(v1-v2));
                if df==0
                    eq=true;
                    if verb
                        fprintf('%s.Text: equal.\n',level);
                    end
                else
                    fprintf('%s.Text: maxdiff %g.\n',level,df);
                    eq=false;
                end
            end
        else
            % Field present in both structures.
            eq=compare(getfield(s1,flds{i}),getfield(s2,flds{i}),...
                       [level,'.',flds{i}],verb);
        end
    end
elseif iscell(s1) && iscell(s2)
    if length(s1)~=length(s2)
        eq=false;
        fprintf('%s: size difference (%d vs. %d).\n',level,length(s1),...
                length(s2));
    else
        for i=1:length(s1)
            eq=compare(s1{i},s2{i},sprintf('%s{%d}',level,i),verb);
        end
    end
elseif ~strcmp(class(s1),class(s2))
    eq=false;
    fprintf('%s: class mismatch.\n',level);
else
    eq=false;
    fprintf('%s: ''%s'' vs. ''%s''.\n',level,s1,s2);
end
