psStub1='camcal7din';
psDir1=fullfile(fileparts(mfilename('fullpath')),'data','weighted','ps',psStub1);
psFile1=fullfile(psDir1,[psStub1,'.psz']);

pmStub='canon7d20mm';
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
        newImName=strrep(imName,fullfile(baseDir,filesep),'images/');
        [p,n,e]=fileparts(newImName);
        if ~exist(fullfile(psDir1,newImName))
            if exist(fullfile(psDir1,p,upper([n,e])))
                newImName=fullfile(p,upper([n,e]));
            end
        end
        prob.images(i).imName=newImName;
    end

    prob.ctrlPts=round(prob.objPts(prob.objPts(:,1)>1000,:));
    prob.ctrlPts(:,5:7)=1e-6;

    % Map ctrl pt ids to 1..4.
    ids=unique([prob.objPts(:,1);prob.ctrlPts(:,1);prob.markPts(:,2)]);
    ids=[ids(ids>1000);ids(ids<=1000)];
    map=sparse(ids,1,1:length(ids));
    prob.ctrlPts(:,1)=full(map(prob.ctrlPts(:,1)));
    prob.objPts(:,1)=full(map(prob.objPts(:,1)));
    prob.markPts(:,2)=full(map(prob.markPts(:,2)));

    % Enforce square pixels.
    prob.job.defCam(4)=prob.job.defCam(5)*prob.job.imSz(1)/prob.job.imSz(2);
    % Enforce no lens distortion.
    prob.job.defCam(6:end)=0;

    % Adjust pixel coordinates.
    prob.markPts(:,3:4)=prob.markPts(:,3:4)+0.5;

    if 0
        % Only keep first 4 images.
        prob.images=prob.images(1:4);
        prob.markPts=prob.markPts(prob.markPts(:,1)<4,:);
    end
end

pmtops(prob,psFile1)

plytoascii(fullfile(psDir1,'unpacked'));

psStub2='camcal7dout';
psDir2=fullfile(fileparts(mfilename('fullpath')),'data','weighted','ps',psStub2);
psFile2=fullfile(psDir2,[psStub2,'.psz']);

pmtops(prob,psFile2)

plytoascii(fullfile(psDir2,'unpacked'));
