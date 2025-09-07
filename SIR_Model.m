clc; close all; clear;

%Membaca file CSV
file = 'Data COVID-19 in Depok City (July 2020).csv';
a = dlmread(file, ';', 1, 0);

%Ekstrak x, y, dan z dari baris 2 - 32
tdat = 1:31;
x = a(1:31, 1);
y = a(1:31, 2);
z = a(1:31, 3);

t0 = 0.0; % nilai waktu t awal (dalam tahun)
tn = 150.0; % nilai waktu t akhir (dalam tahun)
S0 = 1709;% kondisi awal untuk susceptible
I0 = 323; % kondisi awal untuk infectious
R0 = 468; % kondisi awal untuk recovered
beta = 0.000025; % laju infeksi
gamma = 0.05; % laju recovery

% Parameter numerik dan set matriks kosong
ndata = 1000;
S = zeros(ndata,1);
S(1) = S0;
I = zeros(ndata,1);
I(1) = I0;
R = zeros(ndata,1);
R(1) = R0;
t = linspace(t0, tn, ndata);
h = t(2)-t(1);

% Mendefinsikan fungsi
sus = @(t,S,I) -beta.*S.*I;
infect = @(t,S,I) (beta.*S.*I)-(gamma.*I);
recov = @(t,I) (gamma.*I);

% Runge-Kutta Method
for i = 1:ndata-1
  kS1 = sus(t(i),S(i), I(i));
  kI1 = infect(t(i),S(i), I(i));
  kR1 = recov(t(i), I(i));

  kS2 = sus(t(i)+0.5*h, S(i)+0.5*h*kS1, I(i)+0.5*h*kI1);
  kI2 = infect(t(i)+0.5*h, S(i)+0.5*h*kS1, I(i)+0.5*h*kI1);
  kR2 = recov(t(i)+0.5*h, I(i)+0.5*h*kI1);

  kS3 = sus(t(i)+0.5*h, S(i)+0.5*h*kS2, I(i)+0.5*h*kI2);
  kI3 = infect(t(i)+0.5*h, S(i)+0.5*h*kS2, I(i)+0.5*h*kI2);
  kR3 = recov(t(i)+0.5*h, I(i)+0.5*h*kI2);

  kS4 = sus(t(i)+h, S(i)+h*kS3, I(i)+h*kI3);
  kI4 = infect(t(i)+h, S(i)+h*kS3, I(i)+h*kI3);
  kR4 = recov(t(i)+h, I(i)+h*kI3);

  S(i+1)=S(i)+(1/6)*h*(kS1+2*kS2+2*kS3+kS4);
  I(i+1)=I(i)+(1/6)*h*(kI1+2*kI2+2*kI3+kI4);
  R(i+1)=R(i)+(1/6)*h*(kR1+2*kR2+2*kR3+kR4);
end

% Plot hasil
figure(1)
plot(tdat, x, 'x', 'MarkerEdgeColor', 'r', 'LineWidth', 0.5); hold on;
plot(tdat, y, 'x', 'MarkerEdgeColor', 'g', 'LineWidth', 0.5);
plot(tdat, z, 'x', 'MarkerEdgeColor', 'b', 'LineWidth', 0.5);
plot(t, S, 'r', 'LineWidth', 1);
plot(t, I, 'g', 'LineWidth', 1);
plot(t, R, 'b', 'LineWidth', 1);
hold off;
legend('S (data)', 'I (data)', 'R (data)', 'S (model)', 'I (model)', 'R (model)');
xlabel('Waktu (Hari)','FontSize',12)
ylabel('Populasi (Penduduk)','FontSize',12)
title('SIR Model for Covid-19 with 4th Order Runge-Kutta ODE','FontSize',12)
grid on
