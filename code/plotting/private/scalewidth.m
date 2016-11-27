function scalewidth(ax)
% Downscales axes to have same width.


% Scale by c<1 to avoid auto-resize if legend is moved interactively.
c=0.9;

if isscalar(ax)
    p=get(ax,'pos');
else    
    p=cell2mat(get(ax,'pos'));
end
p(:,3)=min(p(:,3))*c;
for i=1:length(ax)
    set(ax(i),'pos',p(i,:));
end
