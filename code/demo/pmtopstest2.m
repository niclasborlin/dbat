psStub1='pmtestin';
psDir1=fullfile(fileparts(mfilename('fullpath')),'data','weighted','ps',psStub1);
psFile1=fullfile(psDir1,[psStub1,'.psz']);

pmStub='weighted-bundle-1mm';
pmFile=fullfile(fileparts(mfilename('fullpath')),'data','weighted','pm',...
                'camcal',pmStub,[pmStub,'-pmexport.txt']);

if ~exist('prob')
    prob=loadpm(pmFile);
    if any(isnan(prob.job.imSz))
        error('No image size');
    end
    baseDir=unique(cellfun(@(x)fileparts(strrep(x,'\',filesep)), ...
                           {prob.images.imName},'uniformoutput', ...
                           false));
    baseDir=baseDir{1};
    for i=1:length(prob.images)
        imName=strrep(prob.images(i).imName,'\',filesep);
        prob.images(i).imName=strrep(imName,...
                                     fullfile(baseDir,filesep),'images/');
    end
    
    % Only keep first 4 images.
    prob.images=prob.images(1:4);
    prob.markPts=prob.markPts(prob.markPts(:,1)<4,:);
end

pmtops(prob,psFile1)

plytoascii(fullfile(psDir1,'unpacked'));

psStub2='pmtestout';
psDir2=fullfile(fileparts(mfilename('fullpath')),'data','weighted','ps',psStub2);
psFile2=fullfile(psDir2,[psStub2,'.psz']);

pmtops(prob,psFile2)

plytoascii(fullfile(psDir2,'unpacked'));
