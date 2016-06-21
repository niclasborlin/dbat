function [prob,err]=loadpm(name)
%LOADPM Load photomodeler export text file into Matlab.
%
%[prob,err]=loadpm(name)
%name - file name to load from.
%prob - struct with fields
%       job      - struct with fields
%                  title   - job title.
%                  tol     - required solution tolerance.
%                  maxIter - maximum number of iterations.
%                  defStd  - vector of default standard deviations 
%                            [photo-coordinates (mm),control points (m),
%                             object points (m),camera position x,y,z (m),
%                             camera angles kappa,phi,omega (degrees)]
%                  defCam  - vector of default calibrated camera
%                            [calibrated focal length (mm), 
%                             principal point xp,yp (mm), 
%                             format size xs,ys (mm),
%                             lens distortion parameters K1,K2,P1,P2]
%                  defCamStd - standard deviation for defCam parameters.
%       cameras  - array of structs with fields
%                  imName   - image name.
%                  outer    - outer orientation parameters 
%                             [x,y,z (m), kappa,phi,omega (degrees)]
%                  outerStd - standard deviation of outer parameters.
%                  outerCov - camera position covariances (xy,xz,yz) (m)
%                  inner    - [focal length, principal point xp,yp, 
%                              format size xs,ys, lens distortion K1,K2,P1,P2]
%                  innerStd - standard deviation of inner parameters.
%       ctrlPts  - matrix of control points [pos x,y,z (m),stdev x,y,z (m)]
%       objPts   - matrix of object points [pos x,y,z (m),stdev x,y,z (m)]
%       markPts  - matrix of marked points [photo#,pt#,pos x,y(px),std x,y(mm)]
%       features - cell array of vectors with points numbers in each feature.
%       featVis  - matrix with rows [i,j] indicating that feature j is 
%                  visible in photo i.

% v1.0  2003-09-17. Niclas Borlin, niclas@cs.umu.se.

err=[];
prob=[];

% Open file.
fid=fopen(name,'rt');
if (fid<0)
	err=['Failed to open file ',name];
	if (nargout<2), error(err); else return; end
end

% Get file size.
if (fseek(fid,0,'eof')~=0)
	err='Failed to get file size';
	if (nargout<2), error(err); else return; end
end
sz=ftell(fid);
% Rewind.
fseek(fid,0,'bof');

h=waitbar(0,'Loading Photomodeler export file');


% Get global information.

% Title
title=fgetl(fid);
% Tolerence, maxIter
s=fgetl(fid); tol=sscanf(s,'%g');
% Default pt stdev.
s=fgetl(fid); defStd=sscanf(s,'%g');
% Default camera params.
s=fgetl(fid); defCam=sscanf(s,'%g');
% Default stdev of camera params.
s=fgetl(fid); defCamStd=sscanf(s,'%g');

job=struct('title',title,'tol',tol(1),'maxIter',tol(2),...
		   'defStd',defStd,'defCam',defCam,'defCamStd',defCamStd);

cameras=[];
while (~feof(fid))
	s=fgetl(fid);
	if (~ishandle(h))
		% Waitbar closed, abort.
		fclose(fid);
		err='Aborted by user';
		if (nargout<2), error(err); else return; end
	else
		waitbarmaybe(sz,ftell(fid)/sz,h);
	end
	% Get photo name.
	[photo,count,err,next]=sscanf(s,'%d');
	if (isempty(photo))
		break;
	end
	imName=s(next:end);
	% Get outer parameters.
	s=fgetl(fid);
	outer=sscanf(s,'%g')';
	outer(1)=[];
	% Get outer stdevs.
	s=fgetl(fid);
	outerStd=sscanf(s,'%g')';
	outerStd(1)=[];
	% Get outer covariances.
	s=fgetl(fid);
	outerCov=sscanf(s,'%g')';
	if (isempty(outerCov))
		outerCov=nan*ones(1,3);
	else
		outerCov(1)=[];
	end
	% Get inner parameters.
	s=fgetl(fid);
	inner=sscanf(s,'%g')';
	inner(1)=[];
	% Get inner stdevs.
	s=fgetl(fid);
	innerStd=sscanf(s,'%g')';
	innerStd(1)=[];
	cameras=[cameras;struct('imName',imName,'outer',outer,...
							'outerStd',outerStd,'outerCov',outerCov,...
							'inner',inner,'innerStd',innerStd)];
end

ctrlPts=[];
while (~feof(fid))
	s=fgetl(fid);
	if (~ishandle(h))
		% Waitbar closed, abort.
		fclose(fid);
		err='Aborted by user';
		if (nargout<2), error(err); else return; end
	else
		waitbarmaybe(sz,ftell(fid)/sz,h);
	end
	% Get control point data [id,x,y,z,sx,sy,sz].
	cp=sscanf(s,'%g')';
	if (isempty(cp))
		break;
	end
	ctrlPts=[ctrlPts;cp];
end

objPts=[];
while (~feof(fid))
	s=fgetl(fid);
	if (~ishandle(h))
		% Waitbar closed, abort.
		fclose(fid);
		err='Aborted by user';
		if (nargout<2), error(err); else return; end
	else
		waitbarmaybe(sz,ftell(fid)/sz,h);
	end
	% Get object point data [id,x,y,z,sx,sy,sz].
	op=sscanf(s,'%g')';
	if (isempty(op))
		break;
	end
	objPts=[objPts;op];
end

markPts=[];
while (~feof(fid))
	s=fgetl(fid);
	if (~ishandle(h))
		% Waitbar closed, abort.
		fclose(fid);
		err='Aborted by user';
		if (nargout<2), error(err); else return; end
	else
		waitbarmaybe(sz,ftell(fid)/sz,h);
	end
	% Get marked point info [ph,mp,x,y,sx,sy].
	mp=sscanf(s,'%g')';
	if (isempty(mp))
		break;
	end
	markPts=[markPts;mp];
end

features=cell(0,0);
while (~feof(fid))
	s=fgetl(fid);
	if (~ishandle(h))
		% Waitbar closed, abort.
		fclose(fid);
		err='Aborted by user';
		if (nargout<2), error(err); else return; end
	else
		waitbarmaybe(sz,ftell(fid)/sz,h);
	end
	% Get pts in each feature.
	ff=sscanf(s,'%d')';
	if (isempty(ff))
		break;
	end
	features{ff(1)}=ff(2+[1:ff(2)]);
end

featVis=[];
while (~feof(fid))
	s=fgetl(fid);
	if (~ishandle(h))
		% Waitbar closed, abort.
		fclose(fid);
		err='Aborted by user';
		if (nargout<2), error(err); else return; end
	else
		waitbarmaybe(sz,ftell(fid)/sz,h);
	end
	% Get photo-feature pair.
	fp=sscanf(s,'%d')';
	if (isempty(fp))
		break;
	end
	featVis=[featVis;fp];
end

prob=struct('job',job,'cameras',cameras,'ctrlPts',ctrlPts,'objPts',objPts,...
			'markPts',markPts,'features',{features},'featVis',featVis);

fclose(fid);
close(h);


function waitbarmaybe(sz,val,h)

persistent call

if (isempty(call) | val==0)
    call=0;
end

call=call+1;
if (rem(call,floor(sqrt(sz/40)))==0)
    waitbar(val,h);
end
    