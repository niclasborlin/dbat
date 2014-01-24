function bundle_result_file(s,e,f)
%BUNDLE_RESULT_FILE Generate result file of bundle run.
%
%   BUNDLE_RESULT_FILE(S,E,F), where S and E are BUNDLE return files and
%   F is a string, writes a text result file to the file F. The text
%   result file contain information about the project, the status of the
%   estimation process, and the quality of the result.
%
%See also: BUNDLE, BUNDLE_COV. 

% $Id$

fid=fopen(f,'wt');
if fid<0
    error('DBAT:bundle_result_file:fileError',...
          ['Failed to open file ''',f,''' for writing.']);
end


% Header info
fprintf(fid,'Damped Bundle Adjustment Toolbox result file\n');
p=repmat(' ',1,3);
fprintf(fid,[p,'Project Name: %s\n'],s.title);

% Info about last bundle.
fprintf(fid,[p,'Information from last bundle\n']);
fprintf(fid,[p,p,'Last Bundle Run: %s\n'],e.dateStamp);
fprintf(fid,[p,p,'DBAT version: %s\n'],e.version);
if e.code==0
    status='OK';
else
    status='fail';
end
fprintf(fid,[p,p,'Status: %s (%d)\n'],status,e.code);
fprintf(fid,[p,p,'Sigma0 (pixels): %g\n'],e.s0);
fprintf(fid,[p,p,'Processing options:\n']);
fprintf(fid,[p,p,p,'Orientation: on\n']);
fprintf(fid,[p,p,p,'Global optimization: on\n']);
offon={'off','on'};
fprintf(fid,[p,p,p,'Calibration: %s\n'],offon{double(any(s.cIO(:)))+1});
fprintf(fid,[p,p,p,'Constraints: off\n']);
fprintf(fid,[p,p,p,'Maximum # of iterations: %d\n'],e.maxIter);
fprintf(fid,[p,p,p,'Convergence tolerance: %g\n'],e.convTol);
fprintf(fid,[p,p,p,'Singular test: %s\n'],offon{double(e.singularTest)+1});
fprintf(fid,[p,p,p,'Chirality veto: %s\n'],offon{double(e.chirality)+1});
fprintf(fid,[p,p,p,'Damping: %s\n'],e.damping.name);

fprintf(fid,[p,p,'Total error:\n']);
fprintf(fid,[p,p,p,'Initial value comment: %s\n'],s.x0desc);
fprintf(fid,[p,p,p,'Number of stages: %d\n'],1);
fprintf(fid,[p,p,p,'Number of iterations: %d\n'],e.usedIters);
fprintf(fid,[p,p,p,'First error: %g\n'],e.res(1));
fprintf(fid,[p,p,p,'Last error: %g\n'],e.res(end));

fprintf(fid,[p,p,'Precisions / Standard Deviations:\n']);
fprintf(fid,[p,p,p,'Photograph Standard Deviations:\n']);


fclose(fid);
