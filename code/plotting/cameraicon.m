function [XYZout,C,body,flash,lens]=cameraicon(v,m)
%CAMERAICON Create camera object for 3D plotting.
%
%[XYZ,C]=cameraicon(v[,m])
%v - vector with camera sizes
%    [w,h,d,l,lr,fr,fd] or
%    [w,h,d] or
%    [s]
%    w  - camera width = -w/2..+w/2. Defaults to 0.11.
%    h  - camera height = -h/2..+h/2. Defaults to 0.08.
%    d  - camera depth = -d..0. Defaults to 0.04;
%    l  - camera lens length = 0..l. Defaults to d/2.
%    lr - lens radius. Defaults to min(w,h)/4.
%    fr - flash radius. Defaults to w/6.
%    fd - flash depth. Defaults to 0.5*d.
%    s  - all sizes are s times their default. Defaults to 1.
%m - should the camera be mirrored in the Z direction? Default: false.
%XYZ,C - array with 3D points such that
%        surf(XYZ(:,:,1),XYZ(:,:,2),XYZ(:,:,3),C) plots the camera.
%        Flash has color==1, lens has color==2, body has color==3.
%
%cameraicon without any return arguments plots the camera.

if (nargin<1), v=1; end
if (nargin<2), m=false; end

if (length(v)==1)
    % Scalar as argument.
    ww=[0.11,0.08,0.04]*v;
else
    % Vector as argument
    ww=v;
end

if (length(ww)<1), w=0.11; else w=ww(1); end
if (length(ww)<2), h=0.08; else h=ww(2); end
if (length(ww)<3), d=0.04; else d=ww(3); end
if (length(ww)<4), l=d/2; else l=w(4); end
if (length(ww)<5), lr=min(w,h)/4; else lr=w(5); end
if (length(ww)<6), fr=w/6; else fr=w(6); end
if (length(ww)<7), fd=d/2; else fd=w(7); end

% Surface separator.
sep=repmat(nan,2,1);

% Unit cube.
boxX=[0,0,1,1,nan,1,1,1,1;
	  0,0,1,1,nan,0,0,0,0];
boxY=[0,1,1,0,nan,1,0,0,1;
	  0,1,1,0,nan,1,0,0,1];
boxZ=[0,0,0,0,nan,1,1,0,0;
	  1,1,1,1,nan,1,1,0,0];

% Scale and position body.
bodyX=(boxX*2-1)*w/2;
bodyY=(boxY*2-1)*h/2;
bodyZ=(boxZ-1)*d;

% Create base cylinder.
[cylX,cylY,cylZ]=cylinder;

% Scale and position.
lensX=cylX*lr;
lensY=cylY*lr;
lensZ=cylZ*l;

% Create lens 'cap'.
lsX=diag([0,1])*lensX;
lsY=diag([0,1])*lensY;
lsZ=lensZ([2,2],:);

% Attach.
lensX=[lensX,sep,lsX];
lensY=[lensY,sep,lsY];
lensZ=[lensZ,sep,lsZ];

% Only keep half the cylinder.
flashX=cylX(:,1:floor(end/2)+1)*fr;
% Put it on top of camera.
flashY=cylY(:,1:floor(end/2)+1)*fr+h/2;
% Make its depth fd*d, and position centered at -d/2.
flashZ=(cylZ(:,1:floor(end/2)+1)*2-1)*fd/2-d/2;

% Create flash 'half-caps'.
fcX=diag([0,1])*flashX;
fcY=diag([0,1])*cylY(:,1:floor(end/2)+1)*fr+h/2;
fcZ1=flashZ([1,1],:);
fcZ2=flashZ([2,2],:);

% Attach half-caps.
flashX=[flashX,sep,fcX,sep,fcX];
flashY=[flashY,sep,fcY,sep,fcY];
flashZ=[flashZ,sep,fcZ1,sep,fcZ2];

if (m)
    bodyZ=-bodyZ;
    lensZ=-lensZ;
    flashZ=-flashZ;
end

% Construct complete 3D array.
XYZ=cat(3,...
		[bodyX,sep,lensX,sep,flashX],...
		[bodyY,sep,lensY,sep,flashY],...
		[bodyZ,sep,lensZ,sep,flashZ]);
    
C=[repmat(3,size(bodyX)+[0,1]),repmat(2,size(lensX)+[0,1]),...
   repmat(1,size(flashX))];

if (nargout<1)
    % No return arguments, visualize camera.
	surf(XYZ(:,:,1),XYZ(:,:,2),XYZ(:,:,3),C);
	colormap(eye(3));
	set(gcf,'renderer','opengl')
	axis equal
	view(3)
	set(gca,'cameraupvector',[0,1,0])
	set(gca,'cameraposition',[-80,30,-20])
	lighting phong
	camlight headlight
	material shiny
	set(gca,'projection','perspective')
	xlabel x
	ylabel y
	zlabel z
else
	XYZout=XYZ;
end

if (nargout>2)
    body=cat(3,bodyX,bodyY,bodyZ);
end

if (nargout>3)
    flash=cat(3,flashX,flashY,flashZ);
end

if (nargout>4)
    lens=cat(3,lensX,lensY,lensZ);
end
