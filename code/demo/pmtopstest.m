psStub='pmtestin';
psDir=fullfile(fileparts(mfilename('fullpath')),'data','weighted','ps',psStub);
psFile=fullfile(psDir,[psStub,'.psz']);

pmStub='weighted-bundle-1cm';
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
end

pmtops(prob,psFile)

plytoascii(fullfile(psDir,'unpacked'));

psStub='pmtestout';
psDir=fullfile(fileparts(mfilename('fullpath')),'data','weighted','ps',psStub);
psFile=fullfile(psDir,[psStub,'.psz']);

pmtops(prob,psFile)

plytoascii(fullfile(psDir,'unpacked'));
