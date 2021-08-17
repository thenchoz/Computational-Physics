function Analyse3(outputFile)

output = load(outputFile);
t = output(:,1);
dt = output(:,2);
x1 = output(:,3);
y1 = output(:,4);
vx1 = output(:,5);
vy1 = output(:,6);
ecin1 = output(:,7);
epot1 = output(:,8);
x2 = output(:,9);
y2 = output(:,10);
vx2 = output(:,11);
vy2 = output(:,12);
ecin2 = output(:,13);
epot2 = output(:,14);
x3 = output(:,15);
y3 = output(:,16);
vx3 = output(:,17);
vy3 = output(:,18);
ecin3 = output(:,19);
epot3 = output(:,20);
ecintot = output(:,21);
epottot = output(:,22);


clear output


figure
plot(x1,y1,'*',x2,y2,'*',x3,y3,'*')
xlabel('x')
ylabel('y')
legend('comète','soleil', 'jupiter')