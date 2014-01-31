function [c,cr,cc,cl,ch]=coverage(s,varargin)
%COVERAGE Percentange of image covered by measurements.
%
%   C=COVERAGE(S), where S is a struct returned by PROB2DBATSTRUCT with N
%   images, returns an N-vector C with the coverage of each image in S. The
%   covered is computed as the area of the convex hull of all measured
%   points in the image, as a percentage of the total area.
%
%   [C,CR]=COVERAGE(S) furthermore returns the rectangular coverage of
%   each image in the N-vector CR.
%
%   [C,CR,CC,CL,CH]=... returns details about each computed coverage. The
%   N-cell array CC contains 2-by-M arrays with the point defining the
%   convex hull as a closed polygon (first point repeated at end). The
%   2-by-N arrays CL=[XMIN;YMIN] and CH=[XMAX,YMAX] contains the minimum and
%   maximum X and Y values, respectively.
%
%   ...=COVERAGE(S,IX), where IX is an N-vector of indices, returns the
%   coverage for the images in IX only. IX='all' corresponds to all images.
%
%   ...=COVERAGE(...,TRUE) treats all considered measurements as being in
%   the same image, i.e. returns one C value, one CR value, etc. Will
%   give unpredictable results if the images have varying sizes.
%
%See also: PROB2DBATSTRUCT, CONVHULL.

% $Id$

ix=1:length(s.imNames);
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

if strcmp(ix,'all'), ix=1:length(s.imNames); end

if union
    % Compute the union of the coverage.
    
    % Assume identical image sizes.
    totArea=prod(s.IO(end-3:end-2,s.cams(ix(1))),1);

    % Extract all points from the requested images.
    i=s.colPos(:,ix);
    pts=s.markPts(:,i(i~=0));

    if isempty(pts)
        cl=nan(2,1);
        ch=nan(2,1);
        cc={pts};
        c=0;
        cr=0;
    else
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
    totArea=prod(s.IO(end-3:end-2,s.cams(ix)),1);

    % Pre-allocate return variables.
    cc=cell(size(ix));
    cl=nan(2,length(ix));
    ch=nan(2,length(ix));
    
    hullArea=nan(size(ix));
    
    for i=1:length(ix)
        % Extract points measured in this image.
        pts=s.markPts(:,s.colPos(s.vis(:,ix(i)),ix(i)));
        
        % Rectangular extents.
        cl(:,i)=min(pts,[],2);
        ch(:,i)=max(pts,[],2);
        
        % Convex hull.
        if nnz(s.vis(:,ix(i)))<3
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
    c=hullArea./totArea;
    cr=prod(ch-cl,1)./totArea;
end
