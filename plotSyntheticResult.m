function plotMDL(skelname,fid,acn,ac1,ac2)
%[k,kk,a0,ssum,MDL,MDL1,MDL2]=load(sprintf('%s.dat',skelname));

aa=load(sprintf('%s.dat',skelname));

figure(1+fid);clf
%p1=plot(aa(:,2),log(aa(:,10)),'+','Color','r')
hold on
p2=plot(aa(:,2),(aa(:,6)),'s','MarkerSize',5,'Color','b')
%p3=plot(aa(:,2),log(aa(:,4)),'o','MarkerSize',3,'Color','b')
title('log(KMCR1)')
xlabel('k')
%legend([p2 p1 p3],'1stage','2stage','residual')
%print('-dpng','testdistKMDL1by2_2stage.png')
plot(acn,ac1,'o','MarkerFaceColor','r','MarkerSize',7)
hold off

figure(2+fid);clf
%p1=plot(aa(:,2),aa(:,11),'+','Color','r')
hold on
p2=plot(aa(:,2),aa(:,7),'s','MarkerSize',5,'Color','b')
hold on
plot(aa(:,2),aa(:,8),'s','MarkerSize',5,'Color','r')
plot(aa(:,2),aa(:,9),'s','MarkerSize',5,'Color','g')
plot(aa(:,2),aa(:,10),'s','MarkerSize',5,'Color','cyan')
plot(aa(:,2),aa(:,11),'s','MarkerSize',5,'Color','k')
plot(acn,ac2(1),'o','MarkerFaceColor','r','MarkerSize',7)
hold off

%p2=plot(aa(:,2),aa(:,12),'s','MarkerSize',5,'Color','b')
title('KMCR2')
xlabel('k')
%legend([p2 p1],'1stage','2stage','Location','northwest')
%print('-dpng','testdistKMDL2by2_2stage.png')

figure(3)
plot(aa(:,6),aa(:,7),'s-','MarkerSize',2,'Color','b')
hold on
plot(aa(:,6),aa(:,8),'s-','MarkerSize',2,'Color','r')
plot(aa(:,6),aa(:,9),'s-','MarkerSize',2,'Color','g')
plot(aa(:,6),aa(:,10),'s-','MarkerSize',2,'Color','cyan')
plot(aa(:,6),aa(:,11),'s-','MarkerSize',2,'Color','k')
plot(ac1,ac2(1),'o','MarkerFaceColor','b','MarkerSize',5)
plot(ac1,ac2(2),'o','MarkerFaceColor','r','MarkerSize',5)
plot(ac1,ac2(3),'o','MarkerFaceColor','g','MarkerSize',5)
plot(ac1,ac2(4),'o','MarkerFaceColor','cyan','MarkerSize',5)
plot(ac1,ac2(5),'o','MarkerFaceColor','k','MarkerSize',5)
hold off

figure(4)
plot(ac2,'o')
%figure(3);clf
%plot(...
%aa(:,9),log
%)

if 1==0
n=3416;
figure(3);clf
for i=1:10
   h=10^-i;
    MDL2=aa(:,2).*log(aa(:,3)/n/h)+n*log(aa(:,4)/n/h+1);
   plot(aa(:,2),MDL2,'+')
   if i==1
    hold on
   end 
end
end
