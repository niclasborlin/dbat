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
fprintf(fid,[p,p,'Sigma0 (pixels): %g\n'],e.s0px);
fprintf(fid,[p,p,'Sigma0 (mm): %g\n'],e.s0mm);
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

corrThreshold=0.95;

if any(s.cIO(:))
    % Camera calibration results.
    fprintf(fid,[p,p,p,'Camera Calibration Standard Deviations:\n']);
    
    % Extract IO covariances and correlations.
    CIO=bundle_cov(s,e,'CIOF');
    [CIOC,ioSigma]=corrmat(CIO,true);
    ioSigma=reshape(ioSigma,size(s.IO,1),[]);

    % Headers and values to print.
    head={'Focal Length','Xp - principal point x','Yp - principal point y',...
          'Fw - format width','Fh - format height',...
          'K1 - radial distortion 1','K2 - radial distortion 2',...
          'K3 - radial distortion 3',...
          'P1 - decentering distortion 1','P2 - decentering distortion 2'};
    names={'Focal','Xp','Yp','Fw','Fh','K1','K2','K3','P1','P2'};
    unit={'mm','mm','mm','mm','mm','mm^(-2)','mm^(-4)','mm^(-6)',...
           'mm^(-2)','mm^(-2)'};
    rows=[3,1:2,11:12,4:6,7:8];    

    S=diag([1,1,-1,1,1,-1,-1,-1,-1,-1]);
    
    for i=1:size(s.cIO,2)
        vals=S*full(s.IO(rows,i));
        sigma=full(ioSigma(rows,i));
        nIO=size(s.IO,1);
        corr=full(CIOC((i-1)*nIO+rows,(i-1)*nIO+rows));
        fprintf(fid,[p,p,p,p,'Camera%d: Unknown\n'],i);
        for j=1:length(head)
            fprintf(fid,[p,p,p,p,p,'%s:\n'],head{j});
            fprintf(fid,[p,p,p,p,p,p,' Value: %g %s\n'],vals(j),unit{j});
            if sigma(j)~=0
                fprintf(fid,[p,p,p,p,p,p,' Deviation: %.1g %s\n'],...
                        sigma(j),unit{j});
            end
            highCorr=find(abs(corr(j,:))>corrThreshold);
            if any(highCorr)
                ss='';
                for k=highCorr(:)'
                    ss=[ss,sprintf(' %s:%.1f%%,',names{k},corr(j,k)*100)];
                end
                ss(end)='.';
                fprintf(fid,[p,p,p,p,p,p,' Correlations over %.1f%%:%s\n'],...
                        corrThreshold*100,ss);
            end
        end    
    end
end

fprintf(fid,[p,p,p,'Photograph Standard Deviations:\n']);

% Get camera station covariances.
CEO=bundle_cov(s,e,'CEOF');
% Compute corresponding correlations and standard deviations.
[CEOC,eoSigma]=corrmat(CEO,true);
eoSigma=reshape(eoSigma,6,[]);

% Headers and values to print.
head={'Omega','Phi','Kappa','Xc','Yc','Zc'};
unit={'deg','deg','deg','','',''};
names=head;
rows=[4:6,1:3];
% Scaling matrix to degrees.
S=diag([180/pi*ones(1,3),ones(1,3)]);

for i=1:size(s.EO,2)
    fprintf(fid,[p,p,p,p,'Photo %d: %s\n'],i,s.imNames{i});
    vals=S*s.EO(rows,i);
    sigma=S*eoSigma(rows,i);
    nEO=size(s.EO,1)-1;
    corr=full(CEOC((i-1)*nEO+rows,(i-1)*nEO+rows));
    for j=1:6
        fprintf(fid,[p,p,p,p,p,'%s:\n'],head{j});
        fprintf(fid,[p,p,p,p,p,p,'Value: %.6f %s\n'],vals(j),unit{j});
        if sigma(j)~=0
            fprintf(fid,[p,p,p,p,p,p,'Deviation: %.1g %s\n'],sigma(j),unit{j});
        end
        highCorr=find(abs(corr(j,:))>corrThreshold);
        if any(highCorr)
            ss='';
            for k=highCorr(:)'
                ss=[ss,sprintf(' %s:%.1f%%,',names{k},corr(j,k)*100)];
            end
            ss(end)='.';
            fprintf(fid,[p,p,p,p,p,p,' Correlations over %.1f%%:%s\n'],...
                    corrThreshold*100,ss);
        end
    end
end

fclose(fid);
