%% SSSP HW 2 ABOUELAZM Youssef - OUALI Ernest
clear; close all; clc

%% Setup

% Simscape File for Ground Truth
load('ssc_output.mat')


fs = 192e3;     % Sampling Frequency [Hz]       
Ts = 1/fs;      % Sampling Period [s]
stop_time = 2;  % Simulation Duration [s]

t = 0:Ts:stop_time; % Time Axis [s]
t = t(1:end-1);

A = 10;             % Signal Amplitude
vin = A * sweeptone(1.99, 0.01, fs, 'SweepFrequencyRange', [500 20000]);    % Input Signal

%% Circuit Parameters
% Piezoelectric MEMS loudspeaker in Free-field conditions

% Transformer Ratio
alpha = 3.7e-4;     % [N / V]
Seff = 2.0e-5;      % [m^2]

% Electrical Domain
Re = 4;             % [Ohms]
Cp = 2.4e-8;        % [F]

% Mechanical Domain
Rm = 9.7e-3;        % [N * s / m]
Mm = 1e-6;          % [kg]
Cm = 2.2e-3;        % [m / N]

% Acoustic Domain
Cbc = 3.6e-13;      % [Pa / m^3]
Ltube1 = 1e2;       % [Pa * s^2 / m^3]
Ltube2 = 1e2;       % [Pa * s^2 / m^3]
Ctube = 6.5e-13;    % [Pa / m^3]
Rac = 5e6;          % [Pa * s / m^3]


%% Removing Ideal Transformers (Mechanical Domain only)
gamma1 = alpha^-1;  % [V / N]
gamma2 = Seff;      % [m^2]

% Resistive Elements
R1 = Re/gamma1^2;   % [Ohms * N^2 / V^2]
R2 = Rm;            % [N * s / m]
R3 = Rac*gamma2^2;  % [Pa * s * m]


% Dynamic Elements
% Capacitors
C1 = Cp * gamma1^2;      % [F * N^2 / V^2]
C2 = Cm;                 % [m / N]
C3 = Cbc / gamma2^2;     % [Pa / m^7]
C4 = Ctube / gamma2^2;   % [Pa / m^7]

% Inductors
L1 = Mm;                 % [kg]
L2 = Ltube1 * gamma2^2;  % [Pa * s^2 * m]
L3 = Ltube2 * gamma2^2;  % [Pa * s^2 * m]

% Source
Fin = alpha .* vin;

%% Adaption Conditions for ports connected to linear elements

Z_R1 = R1;
Z_R2 = R2;
Z_R3 = R3;
Z_C1 = Ts/(2*C1);
Z_C2 = Ts/(2*C2);
Z_C3 = Ts/(2*C3);
Z_C4 = Ts/(2*C4);
Z_L1 = (2*L1)/Ts;
Z_L2 = (2*L2)/Ts;
Z_L3 = (2*L3)/Ts;

% Adaptation conditions for passive elements
Z_2 = Z_R1;
Z_5 = Z_C1;
Z_8 = Z_L1;
Z_11 = Z_R2;
Z_14 = Z_C2;
Z_17 = Z_C3;
Z_20 = Z_L2;
Z_23 = Z_C4;
Z_26 = Z_L3;
Z_27 = Z_R3;

%% Adaptation conditions for adaptors

% Junction I 
Z_25 = Z_26 + Z_27;     % adapted series on port 25


% Junction H: 
Z_24 = Z_25;                            % direct connection with port 25 of junction I
Z_22 = Z_23 * Z_24 / (Z_23 + Z_24);     % adapted parallel on port 22

% Junction G: 
Z_21 = Z_22;                            % direct connection with port 22 of junction H
Z_19 = Z_20 + Z_21;                     % adapted series on port 19

% Junction F: 
Z_18 = Z_19;                            % direct connection with port 19 of junction G
Z_16 = Z_17 * Z_18 / (Z_17 + Z_18);     % adapted parallel on port 16


% Junction E: 
Z_15 = Z_16;                            % direct connection with port 16 of junction F
Z_13 = Z_14 + Z_15;                     % adapted series on port 13

% Junction D: 
Z_12 = Z_13;                            % direct connection with port 13 of junction E
Z_10 = Z_11 + Z_12;                     % adapted series on port 10


% Junction C:
Z_9 = Z_10;                            % direct connection with port 10 of junction D
Z_7 = Z_8 + Z_9;                       % adapted series on port 7

% Junction B: 
Z_6 = Z_7;                             % direct connection with port 7 of junction C
Z_4 = Z_5 * Z_6 / (Z_5 + Z_6);         % adapted parallel on port 4


% Junction A: 
Z_3 = Z_4;                             % direct connection with port 4 of junction B
Z_1 = Z_2 + Z_3;                       % adapted series on port 1

%% Computing Scattering matrices
B = [1 , 1 , 1];
Q = [1 , 1 , 1];

Z_A = diag([Z_1, Z_2, Z_3]);
Z_B = diag([Z_4, Z_5, Z_6]);
Z_C = diag([Z_7, Z_8, Z_9]);
Z_D = diag([Z_10, Z_11, Z_12]);
Z_E = diag([Z_13, Z_14, Z_15]);
Z_F = diag([Z_16, Z_17, Z_18]);
Z_G = diag([Z_19, Z_20, Z_21]);
Z_H = diag([Z_22, Z_23, Z_24]);
Z_I = diag([Z_25, Z_26, Z_27]);


% Series junction: S = 1 - 2 * Z * B' * (B * Z * B')^(-1) * B
S_A = eye(3) - 2 * Z_A * B' * (B * Z_A * B')^(-1) * B;
S_C = eye(3) - 2 * Z_C * B' * (B * Z_C * B')^(-1) * B;
S_D = eye(3) - 2 * Z_D * B' * (B * Z_D * B')^(-1) * B;
S_E = eye(3) - 2 * Z_E * B' * (B * Z_E * B')^(-1) * B;
S_G = eye(3) - 2 * Z_G * B' * (B * Z_G * B')^(-1) * B;
S_I = eye(3) - 2 * Z_I * B' * (B * Z_I * B')^(-1) * B;

% Parallel junction: S = 2 * Q' * (Q * Z^(-1) * Q')^(-1) * Q * Z^(-1) - 1
S_B = 2 * Q' * ( Q * inv(Z_B) * Q' )^(-1) * Q * inv(Z_B) - eye(3);
S_F = 2 * Q' * ( Q * inv(Z_F) * Q' )^(-1) * Q * inv(Z_F) - eye(3);
S_H = 2 * Q' * ( Q * inv(Z_H) * Q' )^(-1) * Q * inv(Z_H) - eye(3);


%% Initialization of Waves

% Incident waves
a1  = 0;   % port 1: incident from source Vin into junction A (Left)
a2  = 0;   % port 2: incident from resistor R1 into junction A (Top)
a3  = 0;   % port 3: incident from junction B into junction A (Right)

a4  = 0;   % port 4: incident from junction A into junction B (Left)
a5  = 0;   % port 5: incident from capacitor C1 into junction B (Top)
a6  = 0;   % port 6: incident from junction C into junction B (Right)

a7  = 0;   % port 7: incident from junction B into junction C (Left)
a8  = 0;   % port 8: incident from inductor L1 into junction C (Top)
a9  = 0;   % port 9: incident from junction D into junction C (Right)

a10 = 0;   % port 10: incident from junction C into junction D (Left)
a11 = 0;   % port 11: incident from resistor R2 into junction D (Top)
a12 = 0;   % port 12: incident from junction E into junction D (Right)

a13 = 0;   % port 13: incident from junction D into junction E (Left)
a14 = 0;   % port 14: incident from capacitor C2 into junction E (Top)
a15 = 0;   % port 15: incident from junction F into junction E (Right)

a16 = 0;   % port 16: incident from junction E into junction F (Left)
a17 = 0;   % port 17: incident from capacitor C3 into junction F (Top)
a18 = 0;   % port 18: incident from junction G into junction F (Right)

a19 = 0;   % port 19: incident from junction F into junction G (Left)
a20 = 0;   % port 20: incident from inductor L2 into junction G (Top)
a21 = 0;   % port 21: incident from junction H into junction G (Right)

a22 = 0;   % port 22: incident from junction G into junction H (Left)
a23 = 0;   % port 23: incident from capacitor C4 into junction H (Top)
a24 = 0;   % port 24: incident from junction I into junction H (Right)

a25 = 0;   % port 25: incident from junction H into junction I (Left)
a26 = 0;   % port 26: incident from inductor L3 into junction I (Top)
a27 = 0;   % port 27: incident from resistor R3 into junction I (Right)

% Reflected waves
b1  = 0;   % port 1: reflected from junction A back to source Vin (Left)
b2  = 0;   % port 2: reflected from junction A into resistor R1 (Top)
b3  = 0;   % port 3: reflected from junction A into junction B (Right)

b4  = 0;   % port 4: reflected from junction B into junction A (Left)
b5  = 0;   % port 5: reflected from junction B into inductor C1 (Top)
b6  = 0;   % port 6: reflected from junction B into junction C (Right)

b7  = 0;   % port 7: reflected from junction C into junction B (Left)
b8  = 0;   % port 8: reflected from junction C into inductor L1 (Top)
b9  = 0;   % port 9: reflected from junction C into junction D (Right)

b10 = 0;   % port 10: reflected from junction D into junction C (Left)
b11 = 0;   % port 11: reflected from junction D into resistor R2 (Top)
b12 = 0;   % port 12: reflected from junction D into junction E (Right)

b13 = 0;   % port 13: reflected from junction E into junction D (Left)
b14 = 0;   % port 14: reflected from junction E into capacitor C2 (Top)
b15 = 0;   % port 15: reflected from junction E into junction F (Right)

b16 = 0;   % port 16: reflected from junction F into junction E (Left)
b17 = 0;   % port 17: reflected from junction F into capacitor C3 (Top)
b18 = 0;   % port 18: reflected from junction F into junction G (Right)

b19 = 0;   % port 19: reflected from junction G into junction F (Left)
b20 = 0;   % port 20: reflected from junction G into inductor L2 (Top)
b21 = 0;   % port 21: reflected from junction G into junction H (Right)

b22 = 0;   % port 22: reflected from junction H into junction G (Left)
b23 = 0;   % port 23: reflected from junction H into capcitor C4 (Top)
b24 = 0;   % port 24: reflected from junction H into junction I (Right)

b25 = 0;   % port 25: reflected from junction I into junction H (Left)
b26 = 0;   % port 26: reflected from junction I into inductor L3 (Top)
b27 = 0;   % port 27: reflected from junction I into resistor R3 (Right)


%% Computational Flow: Forward Scan + Scattering at the Root + Backward Scan

Fout = zeros(size(Fin));


for n = 1 : length(Fin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% UPDATE ACTIVE ELEMENTS %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a5 = b5;        % Capacitor C1
    a8 = -b8;       % Inductor L1
    a14 = b14;      % Capacitor C2
    a17 = b17;      % Capacitor C3
    a20 = -b20;     % Inductor L2
    a23 = b23;      % Capacitor C4
    a26 = -b26;     % Inductor L3


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% FORWARD SCAN %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Junction I
    b25 = S_I(1, :) * [0; a26; a27];  % Scattering relation

    % Junction H
    a24 = b25;                          % Direct transmission from junction I to H
    b22 = S_H(1,:) * [0; a23; a24];   % Scattering relation

    % Junction G
    a21 = b22;                          % Direct transmission from junction I to G
    b19 = S_G(1,:) * [0; a20; a21];   % Scattering relation

    % Junction F
    a18 = b19;                          % Direct transmission from junction G to F
    b16 = S_F(1,:) * [0; a17; a18];   % Scattering relation
    
    % Junction E
    a15 = b16;                          % Direct transmission from junction F to E
    b13 = S_E(1,:) * [0; a14; a15];   % Scattering relation

    % Junction D
    a12 = b13;                          % Direct transmission from junction E to D
    b10 = S_D(1,:) * [0; a11; a12];   % Scattering relation

    % Junction C
    a9 = b10;                           % Direct transmission from junction D to C
    b7 = S_C(1,:) * [0; a8; a9];       % Scattering relation

    % Junction B
    a6 = b7;                            % Direct transmission from junction C to B
    b4 = S_B(1,:) * [0; a5; a6];       % Scattering relation

    % Junction A
    a3 = b4;                            % Direct transmission from junction B to A
    b1 = S_A(1,:) * [0; a2; a3];       % Scattering relation


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% LOCAL ROOT SCATTERING %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    a_root = b1;
    b_root = 2*Fin(n) - a_root;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% BACKWARD SCAN %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Junction A
    a1 = b_root;                        % Direct Transmission from Voltage source to junction A
    b2 = S_A(2,:) * [a1; a2; a3];       % Scattering relation
    b3 = S_A(3,:) * [a1; a2; a3];       % Scattering relation

    % Junction B
    a4 = b3;                            % Direct Transmission from junction A to B
    b5 = S_B(2,:) * [a4; a5; a6];       % Scattering relation
    b6 = S_B(3,:) * [a4; a5; a6];       % Scattering relation

    % Junction C
    a7 = b6;                            % Direct Transmission from junction B to C
    b8 = S_C(2,:) * [a7; a8; a9];       % Scattering relation
    b9 = S_C(3,:) * [a7; a8; a9];       % Scattering relation

    % Junction D
    a10 = b9;                           % Direct Transmission from junction C to D
    b11 = S_D(2,:) * [a10; a11; a12];   % Scattering relation
    b12 = S_D(3,:) * [a10; a11; a12];   % Scattering relation

    % Junction E
    a13 = b12;                          % Direct Transmission from junction D to E
    b14 = S_E(2,:) * [a13; a14; a15];   % Scattering relation
    b15 = S_E(3,:) * [a13; a14; a15];   % Scattering relation

    % Junction F
    a16 = b15;                          % Direct Transmission from junction E to F
    b17 = S_F(2,:) * [a16; a17; a18];   % Scattering relation
    b18 = S_F(3,:) * [a16; a17; a18];   % Scattering relation

    % Junction G
    a19 = b18;                          % Direct Transmission from junction F to G
    b20 = S_G(2,:) * [a19; a20; a21];   % Scattering relation
    b21 = S_G(3,:) * [a19; a20; a21];   % Scattering relation

    % Junction H
    a22 = b21;                          % Direct Transmission from junction G to H
    b23 = S_H(2,:) * [a22; a23; a24];   % Scattering relation
    b24 = S_H(3,:) * [a22; a23; a24];   % Scattering relation

    % Junction I
    a25 = b24;                          % Direct Transmission from junction H to I
    b26 = S_I(2,:) * [a25; a26; a27];   % Scattering relation
    b27 = S_I(3,:) * [a25; a26; a27];   % Scattering relation

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% READ OUTPUT %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fout(n) = (a27 + b27)/2;
end

%% Time-domain Plot

% Computing acoustic pressure
pout = Fout ./ Seff;

% Time Domain Plots, in Pa
figure;
set(gcf, 'Color', 'w');
plot(t, pout, 'b', 'Linewidth', 2);
hold on
plot(gt(1,:), gt(2, :), 'r--', 'Linewidth', 0.5);
grid on;
xlabel('Time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$p_{\mathrm{out}}$ [Pa]','Fontsize',16,'interpreter','latex');
legend('WDF','Ground Truth','Fontsize',16,'interpreter','latex');
title('Output Pressure - Time Domain','Fontsize',18,'interpreter','latex'); 

%% Frequency domain Plot
nfft = 2^20;
res = fs/nfft;
f = (0:nfft/2-1) * res;

ir    = impzest(vin/A, pout);            
ir_gt = impzest(vin/A, gt(2,1:end-1)'); 


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
plot(t, pout - gt(2, 1:end-1)', 'k', 'Linewidth', 2);
grid on;
xlabel('Time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$\mathcal{E}_{\mathrm{out}}$ [Pa]','Fontsize',16,'interpreter','latex');
title('Error Signal','Fontsize',16,'interpreter','latex');

%% Compute Mean Squared Error (MSE)
mse= mean((pout - gt(2, 1:end-1)').^2);
disp('MSE  [Pa^2]= ')
disp(mse)

