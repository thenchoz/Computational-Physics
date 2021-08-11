%<pre><div class='text_to _html'>% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
% 
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour ecraser la valeur d'un parametre par la valeur scannee

%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile
executable = 'Ncorps_OSullivan_Henchoz.e'; % Nom de l'executable
input = 'configuration_corps_8_G1.in'; % Nom du fichier d'entree

nsimul = 8; % Nombre de simulations a faire

precision = logspace(-6,-8,nsimul);
tFin = logspace(2,4, nsimul);

paramstr = 'precision'; % Nom du parametre a scanner
param = precision; % Valeurs du parametre a scanner
paramstr2 = 'tFin'; % Nom du parametre a scanner
param2 = tFin; % Valeurs du parametre a scanner

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, 2 * nsimul);

for ii = 1:nsimul
    filename  = [paramstr , '=', num2str(param (ii))];
    filename2 = [paramstr2, '=', num2str(param2(ii))];
    output{2 * ii - 1} = ['8shaped_', filename , '.out'];
    output{2 * ii    } = ['8shaped_', filename2, '.out'];
    % ! is for Matlab, Octave does not recognize it
    %system(sprintf('%s%s %s %s=%.15g sampling=2 output=%s', repertoire, executable, input, paramstr , param (ii), output{2 * ii - 1}));
    %system(sprintf('%s%s %s %s=%.15g sampling=2 output=%s', repertoire, executable, input, paramstr2, param2(ii), output{2 * ii}));
    disp('Done.')
end
%disp(output)

%% Analyse %%
%%%%%%%%%%%%%

% Parcours des resultats de toutes les simulations

period = zeros(1,nsimul);
period_p = zeros(1,nsimul);
prec1 = zeros(1,nsimul);
prec2 = zeros(1,nsimul);
prec3 = zeros(1,nsimul);

for ii = 1:(2 * nsimul)
  data = load(output{ii});
  
  t = data(:,1);
  if(mod(ii,2) == 0)
    x1 = data(:,3);
    y1 = data(:,4);
    x2 = data(:,9);
    y2 = data(:,10);
  end
  x3 = data(:,15);
  y3 = data(:,16);
  vx3 = data(:,17);
  vy3 = data(:,18);
  
  boucle = 0;
  pos3 = 0;
  size = length(t);
  t_fin = t(size);
  
  for iii = 1:size
    if((x3(iii) < 0) && (y3(iii) < 0) && (vx3(iii) < 0) && (vy3(iii) < 0) && (pos3 == 0))
      ++boucle;
      pos3 = 1;
    elseif((x3(iii) > 0) && (y3(iii) > 0) && (vx3(iii) < 0) && (vy3(iii) < 0) && (pos3 == 1))
      pos3 = 0;
    end
  end
  
  if(mod(ii,2) == 0)
    period_p(ii / 2) = t_fin/(boucle - 1);
  else
    period((ii + 1) / 2) = t_fin/(boucle - 1);
  end
  
  
  if(mod(ii,2) == 0)
    dppp1 = zeros(1,boucle + 1);
    dppp2 = zeros(1,boucle + 1);
    dppp3 = zeros(1,boucle + 1);
    b = 1:boucle + 1;
    boucle = 1;
    
    for iii = 1:size
      ppp1 = sqrt(x1(iii)^2 + y1(iii)^2);
      ppp2 = sqrt(x2(iii)^2 + y2(iii)^2);
      ppp3 = sqrt(x3(iii)^2 + y3(iii)^2);
      
      if(ppp1 < dppp1(boucle) || dppp1(boucle) == 0)
        dppp1(boucle) = ppp1;
      end
      if(ppp2 < dppp2(boucle) || dppp2(boucle) == 0)
        dppp2(boucle) = ppp2;
      end
      if(ppp3 < dppp3(boucle) || dppp3(boucle) == 0)
        dppp3(boucle) = ppp3;
      end
      
      if((x3(iii) < 0) && (y3(iii) < 0) && (vx3(iii) < 0) && (vy3(iii) < 0) && (pos3 == 0))
        ++boucle;
        pos3 = 1;
      elseif((x3(iii) > 0) && (y3(iii) > 0) && (vx3(iii) < 0) && (vy3(iii) < 0) && (pos3 == 1))
        pos3 = 0;
      end
    end
    
    dppp1(boucle + 1) = 0;
    dppp2(boucle + 1) = 0;
    dppp3(boucle + 1) = 0;
    
    prec1(ii / 2) = max(dppp1);
    prec2(ii / 2) = max(dppp2);
    prec3(ii / 2) = max(dppp3);
  end
end


%% Show graphics %%
%%%%%%%%%%%%%%%%%%%

figure
loglog(tFin,period, '*')
xlabel('Temps final')
ylabel('Période')

figure
loglog(precision,period_p, '*')
xlabel('Précision')
ylabel('Période')

figure
loglog(precision,prec1,'-*',precision,prec2,'-*',precision,prec3,'-*')
xlabel('Précision')
ylabel('max des distances min au centre')
legend('C1','C2','C3')


%</div></pre>
