function Analyse_ref_rel_sat(outputFile)

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
%pluton
x11 = output(:,63);
y11 = output(:,64);
vx11 = output(:,65);
vy11 = output(:,66);
ecin11 = output(:,67);
epot11 = output(:,68);
%charon
x12 = output(:,69);
y12 = output(:,70);
vx12 = output(:,71);
vy12 = output(:,72);
ecin12 = output(:,73);
epot12 = output(:,74);
ecintot = output(:,75);
epottot = output(:,76);

clear output

xxJ = x3 - 0;
yyJ = y3 - 0;
distanceJ = sqrt(xxJ.^2 + yyJ.^2);
sinusJ = yyJ ./ distanceJ;
cosinusJ = xxJ ./ distanceJ;

xxT = x6 - x2;
yyT = y6 - y2;
distanceT = sqrt(xxT.^2 + yyT.^2);
sinusT = yyT ./ distanceT;
cosinusT = xxT ./ distanceT;

MPluton = 1.314e22;
MCharon = 1.52e21;
CMx = (x11 * MPluton + x12 * MCharon)/(MPluton + MCharon);
CMy = (y11 * MPluton + y12 * MCharon)/(MPluton + MCharon);
xxP = CMx - x2;
yyP = CMy - y2;
distanceP = sqrt(xxP.^2 + yyP.^2);
sinusP = yyP ./ distanceP;
cosinusP = xxP ./ distanceP;



%FONCTIONNE PAS, CAR TERRE ELLIPSE
%DEUXIEME TENTATIVE
figure %La Terre et la Lune en référentiel relatif
plot(0 ,0,'*',x7.* cosinusT + y7.*sinusT - (x6.* cosinusT + y6.*sinusT) ,y7.*cosinusT - x7.*sinusT -(y6.*cosinusT - x6.*sinusT),'.')
xlabel('x [m]')
ylabel('y [m]')
legend('Terre','Lune')


