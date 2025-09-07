clc; close all; clear;

% Reading CSV File
file = 'Data COVID-19.csv';
a = dlmread(file, ';', 1, 0);

% Extracting Data from CSV for Plotting
tdat = 1:31;
x = a(1:31, 1);
y = a(1:31, 2);
z = a(1:31, 3);

t0 = 0.0;           % Initial t (in days)
tn = 150.0;         % Final t (in days)
S0 = 1709;          % Initial condition for susceptible
I0 = 323;           % Initial condition for infected
R0 = 468;           % Initial condition for recovered
beta = 0.000025;    % Infection rate
gamma = 0.05;       % Recovery rate

% Numerical Setup
ndata = 1000;
S = zeros(ndata,1);
S(1) = S0;
I = zeros(ndata,1);
I(1) = I0;
R = zeros(ndata,1);
R(1) = R0;
t = linspace(t0, tn, ndata);
h = t(2)-t(1);

% Defining the SIR Model Equations using Functions
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

% Plotting the Results
figure(1)
plot(tdat, x, 'x', 'MarkerEdgeColor', 'r', 'LineWidth', 0.5); hold on;
plot(tdat, y, 'x', 'MarkerEdgeColor', 'g', 'LineWidth', 0.5);
plot(tdat, z, 'x', 'MarkerEdgeColor', 'b', 'LineWidth', 0.5);
plot(t, S, 'r', 'LineWidth', 1);
plot(t, I, 'g', 'LineWidth', 1);
plot(t, R, 'b', 'LineWidth', 1);
hold off;
legend('S (data)', 'I (data)', 'R (data)', 'S (model)', 'I (model)', 'R (model)');
xlabel('Time (Days)','FontSize',12)
ylabel('Population (Citizen)','FontSize',12)
title('SIR Model for Covid-19 with 4th Order Runge-Kutta ODE','FontSize',12)
grid on
