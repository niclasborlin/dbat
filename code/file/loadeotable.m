function eo=loadeotable(fName,has)
%LOADEOTABLE Load prior camera positions from text file. 
%
%   EO=LOADEOTABLE(FNAME) loads prior camera positions from the text
%   file FNAME. The data is returned in a struct EO with fields
%       id   - 1-by-M array with id numbers.
%       name - 1-by-M cell array with names.
%       pos  - 3-by-M array of positions.
%       std  - 3-by-M array with prior standard deviations.
%       cov  - 3-by-3-by-M array with prior covariance matrices.
%
%   Each line is expected to contain a comma-separated list of an
%   integer id, a name, X, Y, Z positions, and optional std/covariance
%   information. Blank lines and lines starting with # are ignored.
%
%   EO=LOADEOTABLE(FNAME), where HAS_ID and HAS_NAME are logical,
%   indicates whether the ID and/or NAME are present in the file. At
%   least one of HAS_ID and HAS_NAME must be TRUE.
%
%   The std/covariance information may consist of:
%     - no value - pt is fixed,
%     - 1 value  - sigma_xyz,
%     - 2 values - sigma_xy, sigma_z,
%     - 3 values - sigma_x, sigma_y, sigma_z,
%     - 9 values - full covariance matrix.
%
%NOTE: The functionality is currently identical to that of LOADCPT.

% Call loadcpt to do all the work.
if nargin<2
    eo=loadcpt(fName);
else
    eo=loadcpt(fName,has);
end
