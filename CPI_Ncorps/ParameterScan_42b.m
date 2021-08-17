%<pre><div class='text_to _html'>% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
% 
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour ecraser la valeur d'un parametre par la valeur scannee

%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = ''; % Chemin d'acces au code compile
executable = 'Exercice4.e'; % Nom de l'executable
input = 'configuration_4.2b.in'; % Nom du fichier d'entree

nsimul = 15; % Nombre de simulations a faire

dt = logspace(6,2,nsimul);
sample = 10^6./dt;

xS_fin = 0;
yS_fin = 0;
xC_fin = 4.5442*10^12;
yC_fin = 0;

paramstr = 'dt'; % Nom du parametre a scanner
param = dt; % Valeurs du parametre a scanner

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul);

for ii = 1:nsimul
    filename = [paramstr, '=', num2str(param(ii))];
    output{ii} = ['4.2b', filename, '.out'];
    eval(sprintf('!%s %s %s %s=%.15g sampling=%.15g output=%s', repertoire, executable, input, paramstr, param(ii), sample(ii), output{ii}));
    disp('Done.')
end

%% Analyse %%
%%%%%%%%%%%%%

% Parcours des resultats de toutes les simulations

xS = zeros(1,nsimul);
yS = zeros(1,nsimul);
xC = zeros(1,nsimul);
yC = zeros(1,nsimul);

nstep = zeros(1,nsimul);

Emec = zeros(1,nsimul);

for ii = 1:nsimul
    data = load(output{ii});
    
    xC(ii)=data(end,3);
    yC(ii)=data(end,4);
    xS(ii)=data(end,9);
    yS(ii)=data(end,10);
    
    nstep(ii) = length(data(:,3)) *10; %factor du sampling
    
    Emec(ii) = max(abs(data(end,15)+data(end,16)-data(1,15)-data(1,16))); %A CHANGER MAX DE LA VALEUR ABSOLU
end

%regression lin�aire
index = (dt >= 10000) & (dt <= 1000000);   %# Get the index of the line segment
diff = sqrt((xC-xC_fin).^2 + (yC-yC_fin).^2);
p1 = polyfit(log(dt(index)),log(abs(diff(index))),1);
p2 = polyfit(log(dt(index)),log(abs(Emec(index))),1);
p1(1)
p1(2)
p2(1)
p2(2)
c = linspace(log(dt(end)),log(dt(1)));
%


figure
loglog(dt,abs(xS - xS_fin),'*',dt,abs(xC - xC_fin),'*');hold on; %mettre des valeurs absolus
grid on
xlabel('\Deltat [s]')
ylabel('\Deltax final [m]')
legend('Soleil','Com�te')

figure
loglog(dt,abs(yS - xS_fin),'*',dt,abs(yC - yC_fin),'*')
grid on
xlabel('\Deltat [s]')
ylabel('\Deltay final [m]')
legend('Soleil','Com�te')

figure
plot(xS,yS,'*',xC,yC,'*')
grid on
xlabel('x final [m]')
ylabel('y final [m]')
legend('Soleil','Com�te')

figure % Utiliser pour trouver le dt o� l'erreur est inf�rieur � 10000km
loglog(dt,sqrt((xS-xS_fin).^2 + (yS-yS_fin).^2),'*',dt,sqrt((xC-xC_fin).^2 + (yC-yC_fin).^2),'*');hold on;
loglog(exp(c),exp(p1(1)*c+p1(2)),'r');
grid on
xlabel('\Deltat [s]')
ylabel('\DeltaPosition [m]')
legend('Soleil','Com�te')



figure % Utiliser pour trouver le NSTEP o� l'erreur est inf�rieur � 10000km
loglog(nstep,sqrt((xS-xS_fin).^2 + (yS-yS_fin).^2),'*',nstep,sqrt((xC-xC_fin).^2 + (yC-yC_fin).^2),'*');hold on;
grid on
xlabel('N_{Step}')
ylabel('\DeltaPosition [m]')
legend('Soleil','Com�te')

figure % Convergence de l�nergie m�canique
loglog(dt,Emec,'*');hold on;
loglog(exp(c),exp(p2(1)*c+p2(2)),'r');
grid on
xlabel('\Deltat [s]')
ylabel('\Delta E_{max} [J]')


%</div></pre>
