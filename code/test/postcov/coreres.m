figure(2)
clear
tt={};
for i=1:6
    Z=load(sprintf('result_cores%d.mat',i));
    tt{i}=Z.timeTable(:,[1,2,4]);
end
speedUp=nan(3,5);
for j=1:5
    speedUp(:,j)=mean(tt{1}./tt{j+1},1)';
end

clf
cc=get(gca,'colororder');
cc=cc([1,2,3],:);
hold on
for i=1:size(speedUp,1)
    plot(2:6,speedUp(i,:),'x-','color',cc(i,:))
end
hold off
xlabel('Number of cores')
ylabel('Speedup compared to single core')
legend({'CC','SI','IC'},'location','northwest')
