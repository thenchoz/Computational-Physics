function Analyse12(outputFile)

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
x4 = output(:,21);
y4 = output(:,22);
vx4 = output(:,23);
vy4 = output(:,24);
ecin4 = output(:,25);
epot4 = output(:,26);
x5 = output(:,27);
y5 = output(:,28);
vx5 = output(:,29);
vy5 = output(:,30);
ecin5 = output(:,31);
epot5 = output(:,32);
x6 = output(:,33);
y6 = output(:,34);
vx6 = output(:,35);
vy6 = output(:,36);
ecin6 = output(:,37);
epot6 = output(:,38);
x7 = output(:,39);
y7 = output(:,40);
vx7 = output(:,41);
vy7 = output(:,42);
ecin7 = output(:,43);
epot7 = output(:,44);
x8 = output(:,45);
y8 = output(:,46);
vx8 = output(:,47);
vy8 = output(:,48);
ecin8 = output(:,49);
epot8 = output(:,50);
x9 = output(:,51);
y9 = output(:,52);
vx9 = output(:,53);
vy9 = output(:,54);
ecin9 = output(:,55);
epot9 = output(:,56);
x10 = output(:,57);
y10 = output(:,58);
vx10 = output(:,59);
vy10 = output(:,60);
ecin10 = output(:,61);
epot10 = output(:,62);
x11 = output(:,63);
y11 = output(:,64);
vx11 = output(:,65);
vy11 = output(:,66);
ecin11 = output(:,67);
epot11 = output(:,68);
x12 = output(:,69);
y12 = output(:,70);
vx12 = output(:,71);
vy12 = output(:,72);
ecin12 = output(:,73);
epot12 = output(:,74);
ecintot = output(:,75);
epottot = output(:,76);

clear output

figure
plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,x9,y9,x10,y10,x11,y11,x12,y12)
xlabel('x')
ylabel('y')
legend('Comète','Soleil', 'Jupiter','Io','Callisto','Terre','Lune','Neptune','Mercure','Cérès','Pluton','Charon')

figure
plot(t,vx1,t,vx2,t,vx3,t,vx4,t,vx5,t,vx6,t,vx7,t,vx8,t,vx9,t,vx10,t,vx11,t,vx12)
xlabel('t')
ylabel('vx')
legend('Comète','Soleil', 'Jupiter','Io','Callisto','Terre','Lune','Neptune','Mercure','Cérès','Pluton','Charon')

figure
plot(t,vy1,t,vy2,t,vy3,t,vy4,t,vy5,t,vy6,t,vy7,t,vy8,t,vy9,t,vy10,t,vy11,t,vy12)
xlabel('t')
ylabel('vy')
legend('Comète','Soleil', 'Jupiter','Io','Callisto','Terre','Lune','Neptune','Mercure','Cérès','Pluton','Charon')

figure
plot(t,ecin1,t,ecin2,t,ecin3,t,ecin4,t,ecin5,t,ecin6,t,ecin7,t,ecin8,t,ecin9,t,ecin10,t,ecin11,t,ecin12)
xlabel('t')
ylabel('ecin')
legend('Comète','Soleil', 'Jupiter','Io','Callisto','Terre','Lune','Neptune','Mercure','Cérès','Pluton','Charon')

figure
plot(t,epot1,t,epot2,t,epot3,t,epot4,t,epot5,t,epot6)
xlabel('t')
ylabel('epot')
legend('Comète','Soleil', 'Jupiter','Io','Terre','Lune')

figure
plot(t,dt,'*')
xlabel('t')
ylabel('dt')

figure
plot(t,sqrt((x6 - x2).^2 + (y6 - y2).^2),t,sqrt((x8 - x2).^2 + (y8 - y2).^2),t,sqrt((x9 - x2).^2 + (y9 - y2).^2),t,sqrt((x10 - x2).^2 + (y10 - y2).^2),t,sqrt((x11 - x2).^2 + (y11 - y2).^2),t,sqrt((x12 - x2).^2 + (y12 - y2).^2))
xlabel('t')
ylabel('distance au Soleil')
legend('Terre','Neptune','Mercure','Cérès','Pluton','Charon')

