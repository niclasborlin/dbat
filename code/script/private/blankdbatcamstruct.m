function s=blankdbatcamstruct(withCalInfo)
%BLANKDBATCAMSTRUCT Create a blank camera structure
%
%   S=BLANKDBATCAMSTRUCT returns a DBAT camera structure S with the
%   following fields:
%
%       id         - integer camera id.
%       name       - string with camera name.
%       unit       - string with camera unit.
%       sensor     - comma-separated (width, height) in camera units.
%       image      - comma-separated (width, height) in pixel units.
%       aspectDiff - scalar with 1-pixel aspect ratio.
%       focal      - nominal focal length in camera units.
%       cc         - camera constant in camera units.
%       pp         - comma-separated (x,y) with principal point in camera units.
%       nK         - number of radial coefficients for lens distortion.
%       nP         - number of tangential coefficients for lens distortion.
%       K          - nK-vector with radial lens distortion coefficients.
%       P          - nP-vector with tangential lens distortion coefficients.
%       skew       - scalar with skew.
%       model      - scalar with projection model number.
%       calibrated - logical, true if the camera is calibrated.
%
%   All values are either blank strings or NaN, except nK and nP
%   that are zero.
%
%   S=BLANKDBATCAMSTRUCT(TRUE) furthermore adds a field
%   calibration_info with the following fields:
%
%       calibration_date - string with the calibration timestamp.
%       software         - string with the name of the calibration software.
%       software_version - string with the version of the calibration software.
%       std              - struct with fields cc, pp, aspect, skew, K, P with 
%                          estimated posterior standard deviations.
%       high_corr        - struct with fields
%                          - thres - scalar with correlation threshold
%                          - pairs - string with parameter list, e.g.,
%                                    K1,K2,0.98, with high correlations.
%       calibrated_area  - struct with fields
%                          - type -  string with area type
%                          - points - string with point list x,y;x,y;...
%       cov              - struct with fields
%                          - order    - string with parameter ording
%                          - elements - string with subdiagonal
%                                       non-zero elements i,j,v;i,j,v;...

TODO

%See also: XMLTODBATCAMSTRUCT.

% Create blank camera struct. Struct field names are the same as
% the XML field names.
fieldNames={'id','name','unit','sensor','image','aspectDiff','nK','nP',...
            'focal','model','cc','pp',    'K',       'P',       'skew',...
            'calibrated'};
defaults  ={nan, '',    '',    nan(1,2),nan(1,2),nan,        0,   0,...
            nan,    nan,    nan, nan(1,2),zeros(1,0),zeros(1,0),nan,false}; 

s=cell2struct(defaults,fieldNames,2);

