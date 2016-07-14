file='c:\home\edfdata\pm_exports\1990_original_disk1_cam12.txt';
if (~exist('proj','var'))
    disp('Loading PM export file.');
    proj=loadpm(file);
else
    disp('Using loaded PM export file.');
end

% Extract control point positions...
XYZ=proj.ctrlPts(:,2:4)';
% ...and construct diagonal covariance matrix from standard deviations.
C=diag(reshape(proj.ctrlPts(:,5:7)'.^2,[],1));
cpt=PTGaussian(XYZ,C);

% Extract camera information.
camInner=cat(1,proj.cameras.inner)';
camStd=cat(1,proj.cameras.innerStd)';
cam=PTGaussian(camInner,diag(camStd(:).^2));

% Which camera was used for each image.
camIx=zeros(1,length(proj.cameras));
% Find out how many cameras we have.
uniqueCam=cam(:,1);
camIx(1)=1;
for i=2:size(cam,2)
    probNewCam=zeros(1,size(uniqueCam,2));
    for j=1:size(uniqueCam,2)
        % How different is this cam from previous unique cams?
        probNewCam(j)=probdifferent(cam(:,i),uniqueCam(:,j));
    end
    % If this cam is different from all previous cams, add it to unique
    % cam set.
    if (all(probNewCam>0.9))
        uniqueCam(:,end+1)=cam(:,i);
    end
    [dummy,camIx(i)]=min(probNewCam);
end

% Find the number of pixels in each camera.
px=zeros(2,size(uniqueCam,2));
for i=1:size(uniqueCam,2)
    info=imfinfo(proj.cameras(min(find(camIx==i))).imName);
    px(:,i)=[info.Width;info.Height];
end

uniqueCamPx=[uniqueCam([1,4:5],:);px;uniqueCam([2:3,6:end],:)];

% Number of iterations in simulation.
k=10;

% Photomodeler version.
pmVer=5;

if (1)
    % Do a dry run to test PhotoModeler.
    ctrlptsim(proj,pmVer,0,[],uniqueCamPx,camIx,true);
else
    % Do the actual run.
    camVal=[8.9996 8.2446 6.6 1280 1024 3.9322 3.1916 3.288e-3 -2.158e-5 ...
            0 1.829e-4 1.883e-4];
    [X,C,Xi,Ci]=ctrlptsim(proj,pmVer,k,cpt);
    % Plot the 10-sigma covariance ellipse.
    plot(X,10,'facecolor','b');
    hold on
    plot(C(1:3,:),10,'facecolor','r');
    hold off
    axis equal   
end
