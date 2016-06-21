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

eq=Compare(s1,s2,'s',verb);

eq=ComparePly(s1.document.chunks.chunk.frames.frame.point_cloud,...
              s2.document.chunks.chunk.frames.frame.point_cloud,...
              pszUnpackedDir1,pszUnpackedDir2,...
              's.document.chunks.chunk.frames.frame.point_cloud',verb);


function eq=Compare(s1,s2,level,verb)

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
        elseif ismember(flds{i},numericStrings)
            % Special treatment of numeric string fields.
            if isnumeric(getfield(s1,flds{i}))
                v1=getfield(s1,flds{i});
                v2=getfield(s2,flds{i});
            else
                v1=sscanf(getfield(s1,flds{i}),'%g ');
                v2=sscanf(getfield(s2,flds{i}),'%g ');
            end
            if length(v1)~=length(v2)
                fprintf('%s.%s: length difference (%d vs. %d).\n',level,...
                        flds{i},length(v1),length(v2));
                eq=false;
            elseif isempty(getfield(s1,flds{i}))
                % Silent return.
                eq=true;
            elseif isempty(v1)
                % Numeric processing failed. Do a string comparison instead.
                eq=strcmp(getfield(s1,flds{i}),getfield(s2,flds{i}));
                if eq
                    if verb
                        fprintf('%s.%s: equal.\n',level,flds{i});
                    end
                else
                    fprintf('%s: ''%s'' vs. ''%s''.\n',level,s1,s2);
                end
            else
                df=max(abs(v1-v2));
                if df==0
                    eq=true;
                    if verb
                        fprintf('%s.%s: equal.\n',level,flds{i});
                    end
                else
                    fprintf('%s.%s: maxdiff %g.\n',level,flds{i},df);
                    eq=false;
                end
            end
        else
            % Field present in both structures.
            eq=Compare(getfield(s1,flds{i}),getfield(s2,flds{i}),...
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
            eq=Compare(s1{i},s2{i},sprintf('%s{%d}',level,i),verb);
        end
    end
elseif ~strcmp(class(s1),class(s2))
    eq=false;
    fprintf('%s: class mismatch.\n',level);
elseif ischar(s1)
    eq=strcmp(s1,s2);
    if eq
        if verb
            fprintf('%s: equal.\n',level);
        end
    else
        fprintf('%s: ''%s'' vs. ''%s''.\n',level,s1,s2);
    end
elseif isnumeric(s1) && isnumeric(s2)
    if length(s1)~=length(s2)
        fprintf('%s: length difference (%d vs. %d).\n',level,...
                length(s1),length(s2));
        eq=false;
    else
        df=max(abs(s1-s2));
        if df==0
            eq=true;
            if verb
                fprintf('%s: equal.\n',level);
            end
        else
            fprintf('%s: maxdiff %g.\n',level,df);
            eq=false;
        end
    end    
else
    eq=false;
    fprintf('%s: ''%s'' vs. ''%s''.\n',level,s1,s2);
end

function eq=ComparePly(s1,s2,s1dir,s2dir,level,verb)

if isstruct(s1) && isstruct(s2)
    fld1=fields(s1);
    fld2=fields(s2);
    flds=unique({fld1{:},fld2{:}});
    for i=1:length(flds)
        if xor(isfield(s1,flds{i}),isfield(s2,flds{i}))
            eq=false;
        elseif ismember(flds{i},{'path'})
            % Special treatment of path field.
            f1=getfield(s1,flds{i});
            f2=getfield(s2,flds{i});
            [~,~,d1,~]=ply_read(fullfile(s1dir,f1),'tri');
            [~,~,d2,~]=ply_read(fullfile(s2dir,f2),'tri');
            eq=isequal(d1,d2);
            if eq
                if verb
                    fprintf('%s.path (%s) equal.\n',level,f1);
                end
            else
                eq=Compare(d1,d2,sprintf('(%s)',f1),verb);
            end
        else
            % Field present in both structures.
            eq=ComparePly(getfield(s1,flds{i}),getfield(s2,flds{i}),...
                          s1dir,s2dir,[level,'.',flds{i}],verb);
        end
    end
elseif iscell(s1) && iscell(s2)
    if length(s1)~=length(s2)
        eq=false;
    else
        for i=1:length(s1)
            eq=ComparePly(s1{i},s2{i},s1dir,s2dir,sprintf('%s{%d}',level,i),...
                          verb);
        end
    end
elseif ~strcmp(class(s1),class(s2))
    eq=false;
elseif ischar(s1)
    eq=strcmp(s1,s2);
else
    eq=false;
end
