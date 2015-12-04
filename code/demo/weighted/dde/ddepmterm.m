function rc=ddepmterm(ch)
%DDEPMTERM Terminate DDE channel to Photomodeler.

if ispc
    rc=ddeterm(ch);
else
    rc=1;
end

