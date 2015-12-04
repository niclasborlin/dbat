function [ok,num,str]=ddecmd(ch,cmd,timeOut);

if (nargin<3), timeOut=100; end

if ispc
    % Constants
    NUMERIC=[1,0];
    STRING=[1,1];

    str=ddereq(ch,cmd,STRING,timeOut*1000);
    if (isnumeric(str))
        num={str};
    else
        num=num2cell(sscanf(str,'%g'));
    end
    ok=~(length(num)==0 || num{1}==0);
else
    % Simulate on unix.
    if strcmp(cmd,'GetNextPointID')
        persistent NEXT_PT_ID
        
        if isempty(NEXT_PT_ID)
            NEXT_PT_ID=10000;
        else
            NEXT_PT_ID=NEXT_PT_ID+1;
        end
        ok=true;
        num={NEXT_PT_ID,NEXT_PT_ID};
        str='1';
    else
        ok=true;
        num=1;
        str='1';
    end
end
