function [c,cr,crr,cc,cl,ch,crp]=coverage(s,varargin)
%COVERAGE Percentange of image covered by measurements.
%
%   C=COVERAGE(S), where S is a struct returned by PROB2DBATSTRUCT with N
%   images, returns an N-vector C with the coverage of each image in S. The
%   covered is computed as the area of the convex hull of all measured
%   points in the image, as a fraction of the total area.
%
%   [C,CR,CRR]=COVERAGE(S) furthermore returns the rectangular coverage of
%   each image in the N-vector CR and the radial coverage in the N-vector
%   CRR. The radial coverage is computed as the fraction of the maximum
%   radial distance to the principal point.
%
%   [C,CR,CRR,CC,CL,CH,CRP]=... returns details about each computed
%   coverage. The N-cell array CC contains 2-by-M arrays with the point
%   defining the convex hull as a closed polygon (first point repeated at
%   end). The 2-by-N arrays CL=[XMIN;YMIN] and CH=[XMAX,YMAX] contains the
%   minimum and maximum X and Y values, respectively. The 2-by-2-by-N array
%   CRP contains the principal point in the first column and the point with
%   the maximum radial distance as the second over all 1:N images.
%
%   ...=COVERAGE(S,IX), where IX is an N-vector of indices, returns the
%   coverage for the images in IX only. IX='all' corresponds to all images.
%
%   ...=COVERAGE(...,TRUE) treats all considered measurements as being in
%   the same image, i.e. returns one C value, one CR value, etc. Will
%   give unpredictable results if the images have varying sizes.
%
%See also: PROB2DBATSTRUCT, CONVHULL.


ix=1:length(s.EO.name);
union=false;

for i=1:length(varargin)
    if isnumeric(varargin{i}) || strcmp(varargin{i},'all')
        ix=varargin{i};
    elseif islogical(varargin{i})
        union=varargin{i};
    else
        error('Bad argument');
    end
end

if strcmp(ix,'all'), ix=1:length(s.EO.name); end

if union
    % Compute the union of the coverage.
    
    % Assume identical image sizes.
    totArea=prod(s.IO.sensor.imSize(:,ix(1)));

    % Set up for computation of radial distance.
    
    % Scaling matrix from mm to pixels.
    S=diag([1,-1,1])*diag(1./[s.IO.sensor.pxSize(:,ix(1));1]);
    % Subtract principal point.
    PP=eye(3);
    PP(1:2,3)=s.IO.val(2:3,ix(1));
    % Principal point in pixels.
    pp=euclidean(S*PP*homogeneous(zeros(2,1)));

    % Determine maximum radial distance.
    xx=[1,s.IO.sensor.imSize(1,i)]+0.5*[-1,1];
    yy=[1,s.IO.sensor.imSize(2,i)]+0.5*[-1,1];
    corners=[xx([1,1,2,2]);yy([1,2,2,1])];
    radCorner=sqrt(sum(euclidean(PP\(S\homogeneous(corners))).^2,1));
    maxRad=max(radCorner);
    
    % Extract all points from the requested images.
    i=s.IP.ix(:,ix);
    pts=s.IP.val(:,i(i~=0));

    if isempty(pts)
        cl=nan(2,1);
        ch=nan(2,1);
        cc={pts};
        c=0;
        cr=0;
        crr=0;
        crp=zeros(2,2,0);
    else
        % Radial extent.
        radPts=sqrt(sum(euclidean(PP\(S\homogeneous(pts))).^2,1));
        [ptRadMax,pri]=max(radPts);
        crr=ptRadMax/maxRad;
        crp=[pp,pts(:,pri)];
        
        % Rectangular extents.
        cl=min(pts,[],2);
        ch=max(pts,[],2);
    
        % Convex hull.
        if nnz(i)<3
            % Less that 3 points gives empty convex hull.
            hullArea=0;
            if nargout>2
                % All points are on the border of the convex hull.
                cc={pts};
            end
        else
            % Find convex hull.
            [j,hullArea]=convhull(pts');
            if nargout>2
                % Store hull points if asked for.
                cc={pts(:,j)};
            end
        end
        c=hullArea/totArea;
        cr=prod(ch-cl,1)/totArea;
    end
else
    % Compute individual coverage for each image.
    
    % Image areas.
    totArea=prod(s.IO.sensor.imSize(:,ix),1);

    % Pre-allocate return variables.
    cc=cell(size(ix));
    cl=nan(2,length(ix));
    ch=nan(2,length(ix));
    
    hullArea=nan(size(ix));
    crr=nan(size(ix));

    crp=nan(2,2,length(ix));
    
    for i=1:length(ix)
        % Set up for computation of radial distance for the camera used
        % for image ix(i).
    
        % Scaling matrix from mm to pixels.
        S=diag([1,-1,1])*diag(1./[s.IO.sensor.pxSize(:,ix(i));1]);
        % Subtract principal point.
        PP=eye(3);
        PP(1:2,3)=s.IO.val(2:3,ix(i));
        % Principal point in pixels.
        pp=euclidean(S*PP*homogeneous(zeros(2,1)));

        % Determine maximum radial distance. With the center of the
        % pixels at [1,w] x [1,h], the outermost coordinates in the image
        % becomes [0.5,w+0.5] x [0.5,h+0.5].
        xx=[1,s.IO.sensor.imSize(1,i)]+0.5*[-1,1];
        yy=[1,s.IO.sensor.imSize(2,i)]+0.5*[-1,1];
        corners=[xx([1,1,2,2]);yy([1,2,2,1])];
        radCorner=sqrt(sum(euclidean(PP\(S\homogeneous(corners))).^2,1));
        maxRad=max(radCorner);
        
        % Extract points measured in this image.
        pts=s.IP.val(:,s.IP.ix(s.IP.vis(:,ix(i)),ix(i)));
        
        if ~isempty(pts)
            % Radial extent.
            radPts=sqrt(sum(euclidean(PP\(S\homogeneous(pts))).^2,1));
            [ptRadMax,pri]=max(radPts);
            crr(i)=ptRadMax/maxRad;
            crp(:,:,i)=[pp,pts(:,pri)];
        
            % Rectangular extents.
            cl(:,i)=min(pts,[],2);
            ch(:,i)=max(pts,[],2);
        
            % Convex hull.
            if nnz(s.IP.vis(:,ix(i)))<3
                % Less that 3 points gives empty convex hull.
                hullArea(i)=0;
                if nargout>2
                    % All points are on the border of the convex hull.
                    cc{i}=pts;
                end
            else            
                % Find convex hull.
                [j,hullArea(i)]=convhull(pts');
                if nargout>2
                    % Store hull points if asked for.
                    cc{i}=pts(:,j);
                end
            end
        end
    end
    c=hullArea./totArea;
    cr=prod(ch-cl,1)./totArea;
end
