function ch=ddepminit
%DDEPMINIT Set up a DDE channel to Photomodeler.

if ispc
    ch=ddeinit('PhotoModeler','Data');
    if (ch==0)
        error('Failed to connect to PhotoModeler');
    end
else
    ch=1;
end

