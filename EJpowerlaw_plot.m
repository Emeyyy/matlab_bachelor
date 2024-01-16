x = 0:0.00001:1.2;

y1=x;
y2=x.^5;
y3=x.^15;
y4=x.^50;
y5=x.^Inf;

clf
plot(x,y1,x,y2,x,y3,x,y4,x,y5,'LineWidth',2)
ylim([0 2])
legend('n=1','n=5','n=15','n=50','Bean','Location','northwest')
xlabel('J/J_{C}')
ylabel('E/E_{0}')