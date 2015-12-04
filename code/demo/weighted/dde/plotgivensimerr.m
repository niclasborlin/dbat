function h=plotgivensimerr(fig,year,mag,given,givenId,sim)
%h=plotgivensimerr(fig,year,mag,given,givenId,sim)

figure(fig)
h1=plot(given(1:2,:),mag,'r');
hold on
h2=plot(sim(1:2,:),mag,'b');
x=[given(1,:).mean;sim(1,:).mean;nan(1,size(given,2))];
y=[given(2,:).mean;sim(2,:).mean;nan(1,size(given,2))];
%h3=line(x,y,'color',[0,0.5,0]);
text(given(1,:).mean,given(2,:).mean,int2str(givenId'),'horizontal','center');
legend([h1(1),h2(1)],'Input','Output')
%legend([h1(1),h2(1),h3(1)],'Given','Simulated','Diff')
axis equal
hold off
h=gca;
