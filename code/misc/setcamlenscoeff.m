function s=setcamlenscoeff(s,nK,nP)
%SETCAMLENSCOEFF Set the number of camera lens distortion coefficients.
%
%   S=SETCAMLENSCOEFF(S,NK,NP) sets the number of radial and
%   tangential distortion coefficients in the cameras in S to NK and
%   the NP, respectively. All array sizes are adjusted. As many
%   preset values as possible are retained.
%
%See also: SETCAMMODEL.

oldNK=s.IO.model.nK;
oldNP=s.IO.model.nP;

oldKix=5+(1:s.IO.model.nK);
oldPix=5+s.IO.model.nK+(1:s.IO.model.nP);

newKix=5+(1:nK);
newPix=5+nK+(1:nP);
newN=5+nK+nP;

s.IO.model.nK=nK;
s.IO.model.nP=nP;

if nK<oldNK
    oldKix=oldKix(1:nK);
elseif nK>oldNK
    newKix=newKix(1:oldNK);
end

if nP<oldNP
    oldPix=oldPix(1:nP);
elseif nP>oldNP
    newPix=newPix(1:oldNP);
end

oldIx=[1:5,oldKix,oldPix];
newIx=[1:5,newKix,newPix];

oldVal=s.IO.val;
newVal=nan(newN,size(oldVal,2));
newVal(newIx,:)=oldVal(oldIx,:);

s.IO.val=setfield(s.IO.val,oldIx,newIx,nan);
s.IO.type=setfield(s.IO.type,oldIx,newIx,{''});
s.IO.struct.block=setfield(s.IO.struct.block,oldIx,newIx,nan);

s.prior.IO.val=setfield(s.prior.IO.val,oldIx,newIx,nan);
s.prior.IO.std=setfield(s.prior.IO.std,oldIx,newIx,nan);
s.prior.IO.use=setfield(s.prior.IO.use,oldIx,newIx,false);


function newVal=setfield(oldVal,oldIx,newIx,defVal)
%Adjust size and copy values in oldVal to newVal. Copy rows oldIx from
%old value to rows newIx in new array. Use defVal for any new values.

newVal=repmat(defVal,max(newIx),size(oldVal,2));
newVal(newIx,:)=oldVal(oldIx,:);
