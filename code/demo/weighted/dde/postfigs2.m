% Point id to use for illustration.
id=41;


f90='done_1990_nofw';
f07='done_2007_nofw';

% Load master data files.
Z90=load(f90);
n90=nnz(any(Z90.ptsToRecord));
disp(sprintf('Loaded %s with %d simulations',f90,n90))

Z07=load(f07);
n07=nnz(any(Z07.ptsToRecord));
disp(sprintf('Loaded %s with %d simulations',f07,n07))

% Recorded IDs in each data set
idr90=Z90.idsToRecord;
idr07=Z07.idsToRecord;

% Unpack recorded data
d90=reshape(Z90.ptsToRecord(:,1:n90),3,length(idr90),n90);
d07=reshape(Z07.ptsToRecord(:,1:n07),3,length(idr07),n07);

% Load a posteriori covariances
ap90=load('edfdata/orig/1990_cpids.txt');
ap07=load('edfdata/orig/2007_cpids.txt');

% Checkpoint ids
check90=[...
    5, 11, 21, 31, 41, 49, 57, 63, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74,...
75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91,...
92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 105, 106, 107, 108, 110,...
111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124,...
125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138,...
139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152,...
153, 154, 155, 156, 157, 158, 159, 160, 163, 166, 171, 176, 177, 178,...
179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192,...
193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206,...
207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220,...
221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234,...
235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248,...
249, 250, 251, 252, 253, 254, 255, 256, 260, 263, 266, 271, 273, 274,...
275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288,...
289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302,...
303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316,...
317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330,...
331, 332, 333, 334, 335, 336, 340, 343, 348, 352, 353, 354, 355, 356,...
357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370,...
371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384,...
385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398,...
399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412,...
413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426,...
427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440,...
441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454,...
455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 469, 472, 475, 479,...
484, 488, 489, 492, 493, 494, 497, 502];

check07=[5, 11, 21, 31, 41, 49, 57, 63, 163, 166, 171, 176, 260, 263, 266, ...
         271, 340, 343, 348, 352, 469, 472, 475, 479, 484, 488, 497, 502];

common=intersect(check90,check07);

[dummy,ia,ib]=intersect(Z90.ctrlId,common);
pp90=Z90.pt(:,ia);
[dummy,ia,ib]=intersect(Z07.ctrlId,common);
pp07=Z07.pt(:,ia);

% Samples for this point.
s90=d90(:,Z90.idsToRecord==id,:);
s07=d07(:,Z07.idsToRecord==id,:);

ptsim90=PTGaussian(mean(s90,3),cov(squeeze(s90)'));
ptsim07=PTGaussian(mean(s07,3),cov(squeeze(s07)'));

ptap90=PTGaussian(zeros(3,1),diag(ap90(ap90(:,2)==id,7:9).^2));
ptap07=PTGaussian(zeros(3,1),diag(ap07(ap07(:,2)==id,7:9).^2));

% Specified data for this point.
has90data=ismember(id,Z90.ctrlId);
has07data=ismember(id,Z07.ctrlId);
if has90data
    pt90=Z90.pt(:,Z90.ctrlId==id);
else
    pt90=PTGaussian(nan(3,1));
end
if has07data
    pt07=Z07.pt(:,Z07.ctrlId==id);
else
    pt07=PTGaussian(nan(3,1));
end

ppsim90=PTGaussian(zeros(size(pp90)));
ppsim07=PTGaussian(zeros(size(pp07)));
ppapsim90=PTGaussian(zeros(size(pp90)));
ppapsim07=PTGaussian(zeros(size(pp07)));
for i=1:length(common)
    ii=common(i);
    % Samples for this point.
    s90=d90(:,Z90.idsToRecord==ii,:);
    s07=d07(:,Z07.idsToRecord==ii,:);

    ppsim90(:,i)=PTGaussian(mean(s90,3),cov(squeeze(s90)'));
    ppsim07(:,i)=PTGaussian(mean(s07,3),cov(squeeze(s07)'));

    ppapsim90(:,i)=PTGaussian(zeros(3,1),diag(ap90(ap90(:,2)==ii,7:9).^2));
    ppapsim07(:,i)=PTGaussian(zeros(3,1),diag(ap07(ap07(:,2)==ii,7:9).^2));
end

figure(1)
plot([pt90(1:2),pt07(1:2)],3,{'r-','b-'})
hh1=gca;
axis equal
title(sprintf('Check point %d, 3-sigma for survey uncertainty',id))
legend('1990','2007')

figure(2)
plot([ptsim90(1:2)+ptap90(1:2),ptsim07(1:2)+ptap07(1:2)],3,{'r--','b--'})
hh1(2)=gca;
axis equal
title(sprintf('Check point %d, 3-sigma for photogrammetric uncertainty',id))
legend('1990','2007')

% Give figures 2-4 same axis scaling.
xl=get(hh1,'xlim');
mx=max(cat(2,xl{:}));
mn=min(cat(2,xl{:}));
xl=[mn,mx];

yl=get(hh1,'ylim');
mx=max(cat(2,yl{:}));
mn=min(cat(2,yl{:}));
yl=[mn,mx];

set(hh1,'xlim',xl,'ylim',yl)

print(1,'-depsc','/tmp/chk_pt41src_surv.eps');
print(2,'-depsc','/tmp/chk_pt41src_photo.eps');

figure(3)
hh2=gca;
plot([ptsim90(1:2),ptsim07(1:2)],3,{'r-','b-'})
title(sprintf('Check point %d, 3-sigma due to input control point uncertainty',id))
axis equal

figure(4)
hh2(2)=gca;
plot([ptsim90(1:2).mean+ptap90(1:2),ptsim07(1:2).mean+ptap07(1:2)],3,...
     {'r-.','b-.'});
title(sprintf('Check point %d, 3-sigma due to photogrammetric uncertainty',id))
axis equal

figure(5)
hh2(3)=gca;
plot([ptsim90(1:2)+ptap90(1:2),ptsim07(1:2)+ptap07(1:2)],3,{'b--','r--'});
title(sprintf('Check point %d, 3-sigma for combined error sources',id))
axis equal

% Give figures 2-4 same axis scaling.
xl=get(hh2,'xlim');
mx=max(cat(2,xl{:}));
mn=min(cat(2,xl{:}));
xl=[mn,mx];

yl=get(hh2,'ylim');
mx=max(cat(2,yl{:}));
mn=min(cat(2,yl{:}));
yl=[mn,mx];

set(hh2,'xlim',xl,'ylim',yl)

for i=3:5
    print(i,'-depsc',sprintf('/tmp/chk_pt41src_%c.eps',abs('a')-3+i));
end