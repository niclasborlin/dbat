function genrawfile(year)
%GENRAWFILE
%
%   GENRAWFILE(YEAR) generates a raw text file sim_data_YEAR.txt from stored
%   simulation data.  The sim_data_YEAR.txt file contains a descriptive
%   header followed by the data for each point at each simulation.
%
%   The user is given the choice of loading the data from either the final
%   file 'done_YEAR_nofw.mat' or any of the checkpoint files
%   tmpdata_YEAR_nofw_even.mat and tmpdata_YEAR_nofw_odd.mat, whichever
%   files are available.

% $Id$

if year~=1990 && year~=2007
    error('Year must be 1990 or 2007');
end

% Get setup data.
[idsToSimulate,idsToRecord,ctrlId,pt,dummy,PMid]=load_cpt(year,true);

% Output file.
rawFile=fullfile(pwd,sprintf('sim_data_%d.txt',year));

% Load simulation data.
data=loadsimdata(year,true,'a');

ptsToRecord=data.ptsToRecord;
ptsToSimulate=data.ptsToSimulate;

% Find out how many samples are available, if from a non-complete run.
n=nnz(any(ptsToRecord~=0));

if (n<size(ptsToRecord,2))
    disp(sprintf('Detected %d samples of %d, trimming...',n,size(ptsToRecord,2)));
end

ptsToRecord=ptsToRecord(:,1:n);
ptsToSimulate=ptsToSimulate(:,1:n);

if exist(rawFile,'file')
    disp(['File ''',rawFile,''' exists.']);
    r=input('OK to overwrite (y/n)? ','s');
    if (r(1)~='y')
        disp('Exiting without writing anything.')
        return;
    end
end

disp(' ');
fid=fopen(rawFile,'wt+');
if fid<0
    error('Failed to open raw data file for writing.');
else
    disp(['Writing to file ',rawFile,'...']);
    fprintf(fid,'# iteration, id, flag (1=simulated, 0=calculated), X, Y, Z\n');
end

for loop=1:n
    if rem(loop-1,floor(sqrt(n)))==0
        disp(sprintf('Writing loop %d of %d...',loop,n))
    end
    % Print simulation input data to file.
    data=reshape(ptsToSimulate(:,loop),3,[]);
    id=idsToSimulate;
    for i=1:length(id)
        fprintf(fid,'%4d, %3d, 1, %.4f, %.4f, %.4f\n',loop,id(i),data(:,i));
    end
    
    % Write simulation output data to file.
    data=reshape(ptsToRecord(:,loop),3,[]);
    id=idsToRecord;
    for i=1:length(id)
        fprintf(fid,'%4d, %3d, 0, %.4f, %.4f, %.4f\n',loop,id(i),data(:,i));
    end
end
fclose(fid);
disp('Done.')
