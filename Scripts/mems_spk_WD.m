clear; close all; clc

%% Simscape File for Ground Truth

load('ssc_output.mat')

%% Sampling Frequency
fs = 192e3;

%% Sampling Period
Ts = 1/fs;

%% Simulation Duration
stop_time = 2;  % [seconds]

%% Input Signal

% Time Axis
t = 0:Ts:stop_time;
t = t(1:end-1);

% Signal Amplitude
A = 10;
% vin = A * sweeptone(1.99, 0.01, fs, 'SweepFrequencyRange', [500 20000]);
load('vin.mat');

%% Circuit Parameters
% Piezoelectric MEMS loudspeaker in Free-field conditions

% Transformer Ratio
alpha = 3.7e-4;
Seff = 2.0e-5;

% Electrical Domain
Re = 4;
Cp = 2.4e-8;

% Mechanical Domain
Rm = 9.7e-3;
Mm = 1e-6;
Cm = 2.2e-3;

% Acoustic Domain
Cbc = 3.6e-13; 
Ltube1 = 1e2;
Ltube2 = 1e2;
Ctube = 6.5e-13;
Rac = 5e6;

%% Removing Ideal Transformers (Mechanical Domain only)

gamma1 = alpha^-1;
gamma2 = Seff;

% Resistive Elements
R1 = Re/gamma1^2;
R2 = Rm;
R3 = Rac*gamma2^2;


% Dynamic Elements
% Capacitors
C1 = Cp*gamma1^2;
C2 = Cm;
C3 = Cbc/gamma2^2;
C4 = Ctube/gamma2^2;


% Inductors
L1 = Mm;
L2 = Ltube1*gamma2^2;
L3 = Ltube2*gamma2^2;

% Source
Fin = alpha * vin;

%% Setting of Free Parameters (Adaptation Conditions)


%% Computing Scattering Matrices


%% Initialization of Waves


%% Initialization of Output Signals

Fout = zeros(1, length(t));

%% Simulation Algorithm

for n = 1 : length(Fin)

    % Forward Scan
    
    % Local Root Scattering

    % Backward Scan

    % Read Output

end

%% Output Plots

% Computing acoustic pressure
pout = Fout ./ Seff;

% Time Domain Plots
figure
set(gcf, 'Color', 'w');
plot(t, pout, 'b', 'Linewidth', 2);
hold on
plot(gt(1, :), gt(2, :), 'r--', 'Linewidth', 2);
grid on;
xlabel('Time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$p_{\mathrm{out}}$ [Pa]','Fontsize',16,'interpreter','latex');
legend('WDF','SSC','Fontsize',16,'interpreter','latex');
title('Output Pressure - Time Domain','Fontsize',18,'interpreter','latex');


% Frequency domain Plots
nfft = 2^20;
res = fs/nfft;
f = (0:nfft/2-1) * res;

ir = impzest(vin/A, pout');
ir_gt = impzest(vin/A, gt(2, 1:end-1)');

tf = fft(ir, nfft);
tf_gt = fft(ir_gt, nfft);

abs_tf = abs(tf(1:nfft/2));
abs_tf_gt = abs(tf_gt(1:nfft/2));

figure
set(gcf, 'Color', 'w');
semilogx(f, 20*log10(abs_tf/2e-5), 'b', 'Linewidth', 2);
hold on
semilogx(f, 20*log10(abs_tf_gt/2e-5), 'r--', 'Linewidth', 2);
grid on;
xlim([500, 20000])
xlabel('Frequency [Hz]','Fontsize',16,'interpreter','latex');
ylabel('$\mathrm{SPL}\,[\mathrm{dB}_\mathrm{SPL}]$','Fontsize',16,'interpreter','latex');
legend('WDF','SSC','Fontsize',16,'interpreter','latex');
title('Output Sound Pressure Level - Frequency Domain','Fontsize',16,'interpreter','latex');

%% Error Plots

figure
set(gcf, 'Color', 'w');
hold on;
plot(t, pout - gt(2, 1:end-1), 'k', 'Linewidth', 2);
grid on;
xlabel('Time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$\mathcal{E}_{\mathrm{out}}$ [Pa]','Fontsize',16,'interpreter','latex');
title('Error Signal','Fontsize',16,'interpreter','latex');

%% Compute Mean Squared Error (MSE)

mse = mean((pout' - gt(2, 1:end-1)).^2);
disp('MSE = ')
disp(mse)
