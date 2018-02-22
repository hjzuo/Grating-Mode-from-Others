figure(1)

plot(a(1,:),a(2:4,:),'linewidth',2);

legend('610nm','625nm','640nm');

ylabel('GVD(ps*km^-^1nm^-^1)')

xlabel('Wavelength(\mum)');

title('580nm Film Thickness')

axis([1.40,1.45,-30,20])


xx = xlim;

yy = ylim;

line([xx(1),xx(2)],[0,0],'color','c','linewidth',2);

line([xx(1),xx(2)],[7,7],'color','b','linewidth',2);

line([1.42,1.42],[yy(1),yy(2)],'color','k','linewidth',2)

box on