function Analyse3_8shape(outputFile)

%% Lecture inputFile %%
%%%%%%%%%%%%%%%%%%%%%%%

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

%% Analyse rotation %%
%%%%%%%%%%%%%%%%%%%%%%

boucle = -1;
pos3 = 0;
size = length(t);
t_fin = t(size);

for ii = 1:size
  if((x3(ii) < 0) && (y3(ii) < 0) && (vx3(ii) < 0) && (vy3(ii) < 0) && (pos3 == 0))
    ++boucle;
    pos3 = 1;
  elseif((x3(ii) > 0) && (y3(ii) > 0) && (vx3(ii) < 0) && (vy3(ii) < 0) && (pos3 == 1))
    pos3 = 0;
  end
end

disp(boucle)
disp(" passages, T approx ")
disp(t_fin/boucle)

dppp1 = zeros(1,boucle + 1);
dppp2 = zeros(1,boucle + 1);
dppp3 = zeros(1,boucle + 1);
b = 1:boucle + 1;
boucle = 1;

for ii = 1:size
  ppp1 = sqrt(x1(ii)^2 + y1(ii)^2);
  ppp2 = sqrt(x2(ii)^2 + y2(ii)^2);
  ppp3 = sqrt(x3(ii)^2 + y3(ii)^2);
  
  if(ppp1 < dppp1(boucle) || dppp1(boucle) == 0)
    dppp1(boucle) = ppp1;
  end
  if(ppp2 < dppp2(boucle) || dppp2(boucle) == 0)
    dppp2(boucle) = ppp2;
  end
  if(ppp3 < dppp3(boucle) || dppp3(boucle) == 0)
    dppp3(boucle) = ppp3;
  end
  
  if((x3(ii) < 0) && (y3(ii) < 0) && (vx3(ii) < 0) && (vy3(ii) < 0) && (pos3 == 0))
    ++boucle;
    pos3 = 1;
  elseif((x3(ii) > 0) && (y3(ii) > 0) && (vx3(ii) < 0) && (vy3(ii) < 0) && (pos3 == 1))
    pos3 = 0;
  end
end

%% Show graphics %%
%%%%%%%%%%%%%%%%%%%

figure
plot(x1,y1,'-*',x2,y2,'-*',x3,y3,'-*')
xlabel('x')
ylabel('y')
legend('C1','C2','C3')

figure
plot(b,dppp1,'-*',b,dppp2,'-*',b,dppp3,'-*')
xlabel('nb boucle')
ylabel('Distance min au centre')
legend('C1','C2','C3')