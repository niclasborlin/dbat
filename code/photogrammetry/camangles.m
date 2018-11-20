function a=camangles(s,msg)
%CAMANGLES Camera ray angles for a project.
%
%   A=CAMANGLES(S), where S is a struct returned by PROB2DBATSTRUCT
%   with N cameras, returns an N-vector A with the maximum angle in
%   radians between rays for each camera. The maximum angle is the
%   angle closest to being orthogonal between pairs of rays for each
%   camera.
%
%   A zero angle is returned for single-ray cameras. NaN's are
%   returned for cameras without rays.
%
%   A=CAMANGLES(S,MSG), uses MSG as the message for a delayed
%   waitbar. The waitbar is presented if the angle computation takes
%   more than one second.

if nargin<2, msg=''; end

a=zeros(1,size(s.EO.val,2));

% Delayed progress dialog.
start=clock;
lapTime=start;
h=[];

for i=1:length(a)
    % Camera center.
    c=s.EO.val(1:3,i);
    % Object points.
    pp=s.OP.val(:,s.IP.vis(:,i));
    % How many rays?
    switch size(pp,2)
      case 0
        % Do nothing.
      case 1
        a(i)=0;
      otherwise
        % Direction vectors.
        d=repmat(c,1,size(pp,2))-pp;
        % Normalize
        dn=d./repmat(sqrt(sum(d.^2,1)),3,1);
        % Compute all inner products. Guard for round-off errors.
        ip=max(min(dn'*dn,1),-1);
        % Angle is acos of inner product of normalized vectors.
        a(i)=max(acos(abs(ip(:))));
    end
    
    if ~isempty(msg)
        % Use waitbar only if we have a message.
        if isempty(h) && etime(clock,start)>1 && i~=length(a)
            % Only create dialog if execution takes more than 1s and this
            % iteration is not the last.
            h=waitbar(i/length(a),msg);
            lapTime=clock;
        elseif etime(clock,lapTime)>1
            % Update dialog every 1 s.
            if ishandle(h) % Guard against window close.
                waitbar(i/length(a),h);
            end
            lapTime=clock;
        end
    end
end

if ishandle(h), close(h), end
