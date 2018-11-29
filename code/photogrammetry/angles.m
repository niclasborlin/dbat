function a=angles(s,msg)
%ANGLES Object point ray angles for a project.
%
%   A=ANGLES(S), where S is a struct returned by PROB2DBATSTRUCT with N
%   object points, returns an N-vector A with the maximum angle in radians
%   between rays for each object point. The maximum angle is the angle
%   closest to being orthogonal between pairs of rays for each
%   object point.
%
%   A zero angle is returned for single-ray points. NaN's are
%   returned for points without rays.
%
%   A=ANGLES(S,MSG), uses MSG as the message for a delayed
%   waitbar. The waitbar is presented if the angle computation takes
%   more than one second.

if nargin<2, msg=''; end

a=nan(size(s.OP.id));

% Delayed progress dialog.
start=clock;
lapTime=start;
h=[];

for i=1:length(s.OP.id)
    % Point
    p=s.OP.val(:,i);
    % Camera centers.
    cc=s.EO.val(1:3,s.IP.vis(i,:));
    % How many rays?
    switch size(cc,2)
      case 0
        % Do nothing.
      case 1
        a(i)=0;
      otherwise
        % Direction vectors.
        d=repmat(p,1,size(cc,2))-cc;
        % Normalize
        dn=d./repmat(sqrt(sum(d.^2,1)),3,1);
        % Compute all inner products. Guard for round-off errors.
        ip=max(min(dn'*dn,1),-1);
        % Angle is acos of inner product of normalized vectors.
        a(i)=max(acos(abs(ip(:))));
    end
    
    if ~isempty(msg)
        % Use waitbar only if we have a message.
        if isempty(h) && etime(clock,start)>1 && i~=length(s.OP.id)
            % Only create dialog if execution takes more than 1s and this
            % iteration is not the last.
            h=waitbar(i/length(s.OP.id),msg);
            lapTime=clock;
        elseif etime(clock,lapTime)>1
            % Update dialog every 1 s.
            if ishandle(h) % Guard against window close.
                waitbar(i/length(s.OP.id),h);
                pause(0.01);
            end
            lapTime=clock;
        end
    end
end

if ishandle(h), close(h), end
