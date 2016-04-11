psStub='pmtest';
psDir=fullfile(fileparts(mfilename('fullpath')),'data','weighted','ps',psStub);
psFile=fullfile(psDir,[psStub,'.psz']);

pmStub='weighted-bundle-1cm';
pmFile=fullfile(fileparts(mfilename('fullpath')),'data','weighted','pm',...
                'camcal',pmStub,[pmStub,'-pmexport.txt']);

if ~exist('prob')
    prob=loadpm(pmFile);
    baseDir=unique(cellfun(@(x)fileparts(strrep(x,'\',filesep)), ...
                           {prob.images.imName},'uniformoutput',false));
    for i=1:length(prob.images)
        prob.images(i).imName=strrep(prob.images(i).imName,...
                                     fullfile(baseDir,filesep),'images/');
    end
end

pmtops(prob,psFile)
