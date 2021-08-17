function ex2_HCL_interferometer(n)
% ex2_HCL_interferometer : switch between the differents parts of the second ex/report 1.
%
% Arguments:
%       input interger: number of exercice;
%
% Returns : nothing.

data = importdata('Hcl_interferogram.dat');
Delta_D = data(:,1);
With_S  = data(:,2);
Without = data(:,3);

switch n
    case 1
        ex2_1(Delta_D, With_S, Without)           
    case 2
        ex2_2(Delta_D, With_S, Without)
    case 3
        ex2_3(Delta_D, With_S, Without)
    case 4
        ex2_4(Delta_D, With_S, Without)
    case 5
        ex2_5()
end
end

function [N, Delta_x, With_S, Without]=analysed(Delta_D, With_S, Without)
% function that avoid repetition

min_D = min(Delta_D);
max_D = max(Delta_D);
a = max_D - min_D;
N = length(Delta_D);
Delta_x = a / N;

if(min_D == -max_D)
    With_S  = ifftshift(With_S );
    Without = ifftshift(Without);
end
end


function ex2_1(Delta_D, With_S, Without)
% ex2_1 : plot data from Hcl_interferogram.dat

figure
plot(Delta_D, With_S, '.', Delta_D, Without, '.');
legend('With sample', 'Without sample')
xlabel('$\Delta$D [cm]')
ylabel('Intensity')
end

function ex2_2(Delta_D, With_S, Without)
% ex2_2 : convert data with DFT, then plot

[N, Delta_x, With_S, Without] = analysed(Delta_D, With_S, Without);

% con_With_S  = abs(fftshift(mydft(With_S )) .* Delta_x).^2;
% con_Without = abs(fftshift(mydft(Without)) .* Delta_x).^2;
% fft instead of mydft for speed improve
con_With_S  = abs(fftshift(fft(With_S )) .* Delta_x).^2;
con_Without = abs(fftshift(fft(Without)) .* Delta_x).^2;
frequency   = 1/(N * Delta_x) * (-N/2:1:N/2-1);

figure
plot(frequency, con_With_S, frequency, con_Without);
%xlim([2500 3200])
legend('With sample', 'Without sample')
xlabel('$\nu$ [cm$^{-1}$]')
ylabel('$\vert$Intensity$\vert^2$')
end

function ex2_3(Delta_D, With_S, Without)
% ex2_3 : convert data with DFT, then substract & reduce analysed region

[N, Delta_x, With_S, Without] = analysed(Delta_D, With_S, Without);

% substract  = abs(fftshift(mydft(Without)) .* Delta_x - fftshift(mydft(With_S )) .* Delta_x).^2;
substract  = abs(fftshift(fft(Without)) .* Delta_x - fftshift(fft(With_S )) .* Delta_x).^2;
frequency  = 1/(N * Delta_x) * (-N/2:1:N/2-1);

figure
plot(frequency, substract);
xlim([2500 3200])
xlabel('$\nu$ [cm$^{-1}$]')
ylabel('$\vert$Intensity$\vert^2$')
end

function ex2_4(Delta_D, With_S, Without)
% ex2_4 : analysed the J (0,1) peaks (around 2860 & 2905 cm-1) (deduce from ex_3)

[N, Delta_x, With_S, Without] = analysed(Delta_D, With_S, Without);

% substract  = abs(fftshift(mydft(Without)) .* Delta_x - fftshift(mydft(With_S )) .* Delta_x).^2;
substract  = abs(fftshift(fft(Without)) .* Delta_x - fftshift(fft(With_S )) .* Delta_x).^2;
frequency  = 1/(N * Delta_x) * (-N/2:1:N/2-1);

x_i_p1 = find(abs(frequency - 2861) < 0.4); % begining of first pic
x_m_p1 = find(abs(frequency - 2863) < 0.4); % switch first pic
x_f_p1 = find(abs(frequency - 2866) < 0.4); % end of first pic
x_i_p2 = find(abs(frequency - 2902) < 0.4); % begining of second pic
x_m_p2 = find(abs(frequency - 2904) < 0.4); % switch first pic
x_f_p2 = find(abs(frequency - 2907) < 0.4); % end of second pic
% Value optain from graphic, 0.4 for presicion because frequency define
% with Delta_x, a not an integer

size_p11 = 0;
size_p12 = 0;
size_p21 = 0;
size_p22 = 0;
for i = x_i_p1+1:x_m_p1
    size_p11 = size_p11 + (substract(i-1)+substract(i))/2; % trapez integration, manualy done, delta_x = 1
end
for i = x_m_p1+1:x_f_p1
    size_p12 = size_p12 + (substract(i-1)+substract(i))/2; % trapez integration, manualy done, delta_x = 1
end
for i = x_i_p2+1:x_m_p2
    size_p21 = size_p21 + (substract(i-1)+substract(i))/2; % trapez integration, manualy done, delta_x = 1
end
for i = x_m_p2+1:x_f_p2
    size_p22 = size_p22 + (substract(i-1)+substract(i))/2; % trapez integration, manualy done, delta_x = 1
end

% substract background
% background = (substract(x_i_p1) + substract(x_f_p1) + substract(x_i_p2) + substract(x_f_p2)) / 4; % approx background
background = 0;
for i = x_f_p1:x_i_p2
    background = background + substract(i);
end
background = background / (x_i_p2 - x_f_p1 + 1);
size_p11 = size_p11 - background * (x_m_p1 - x_i_p1);
size_p12 = size_p12 - background * (x_f_p1 - x_m_p1);
size_p21 = size_p21 - background * (x_m_p2 - x_i_p2);
size_p22 = size_p22 - background * (x_f_p2 - x_m_p2);

disp('pic 1 : ');
disp('   pic 1.1 : ');
disp(size_p11);
disp('   pic 1.2 : ');
disp(size_p12);
disp('ratio (p1.1/(p1.1+p1.2)) :');
disp(size_p11/(size_p11+size_p12));
disp('pic 2 : ');
disp('   pic 2.1 : ');
disp(size_p21);
disp('   pic 2.2 : ');
disp(size_p22);
disp('ratio (p2.1/(p2.1+p2.2)) :');
disp(size_p21/(size_p21+size_p22));

figure
plot(frequency, substract);
xlim([2858 2912])
xlabel('$\nu$ [cm$^{-1}$]')
ylabel('$\vert$Intensity$\vert^2$')
end

function ex2_5()
% ex2_5 : finding bond lenth

% formula : delta nu = hbar J (J + 1) / (2 pi c mu d^2), mu = MH MCl / (MH + MCl)
% => d = sqrt(hbar J (J + 1) / (2 pi c mu deltanu))

h = 1.054571726 * 10^(-34); % hbar Js
c = 29979245800; % light speed cm/s
J = 1; % Jbar, cause we take 0<->1

amu = 1.660538921 * 10^(-27); % 1 a.m.u. to Kg
MH = 1.0078 * amu; % MH Kg * convert
MCl_35 = 34.969 * amu; % MCl_35 Kg * convert
MCl_37 = 36.966 * amu; % MCl_37 Kg * convert
mu_35 = MH * MCl_35 / (MH + MCl_35);
mu_37 = MH * MCl_37 / (MH + MCl_37);

Delta_nu = 2903-2862; % According to ex2_4 -> pic to pic, small & big have the same diff between them
A = 10^(-10); % 1 °A = 10^-10 m

d_35 = sqrt(h * J * (J + 1) / (2 * pi * c * mu_35 * Delta_nu));
d_37 = sqrt(h * J * (J + 1) / (2 * pi * c * mu_37 * Delta_nu));
disp('d (Cl 35) :');
disp(d_35 / A);
disp('d (Cl 37) :');
disp(d_37 / A);
end
