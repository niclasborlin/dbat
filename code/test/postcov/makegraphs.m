Z1=load('result_pixel');
Z2=load('result_pixel_self');

ix=[1,5,3,2,6:11]

spBix=10;
spLBix=11;
spLCix=12;
spBix=11;

Z1.timeTable(ix,[10:12])

methods={'CC','SI','IC','LD'};
mIx=[1,2,4,5];

figure(1)

loglog(Z1.timeTable(ix,spBix),Z1.timeTable(ix,mIx),'x-')
cc=get(gca,'colororder')
hold on
for i=1:length(mIx)
    loglog(Z2.timeTable(ix,spBix),Z2.timeTable(ix,mIx(i)),'o--','color',cc(i,:))
end
hold off
legend(methods,'location','northwest')

xlabel('Density')
ylabel('Time')

set(gca,'xdir','reverse')
ylim([-inf,400])

