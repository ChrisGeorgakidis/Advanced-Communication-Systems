%*************** Advanced Communication Systems*****************%
%                       CE542, Fall 2018                        %
%                       ECE, UTH, Greece                        %
% File: project_4.m                                            %
% Authors: Christos Georgakidis (1964)                          %
%***************************************************************%

% clear data from previous run %
clear
close all
clc

% System Specifications %
OVERSAMPLE    = 1;
M             = 2;              % Modulation type %
k             = log2(M);
SNRdB         = -50:5:50;           % SNR in dB %
N             = k*10000;        % number of bits, it is multiplied by log2(M) in%
                                % order to ensure that it will be a multiple of it %
Eb            = 1;              % Energy / bit %
SNR_linear    = 10.^(SNRdB/10); % convert SNR from db to linear: SNR(dB) = 10log10(SNR(linear)) %
nOfIterations = 10;             % num of iterations for Monte Carlo Simulation %
T = 5 * k; % The time the channels remains the same %
taps = 5;   % Number of taps %

fprintf('----------------------------------\n');
fprintf('| Discrete LTI Channel with ISI  |\n');
fprintf('| Author: Georgakidis Christos   |\n');
fprintf('| Date: 20/12/2018               |\n');
fprintf('----------------------------------\n\n');

fprintf('---------------------------------------\n');
fprintf('Discrete LTI Channel with ISI equaliser\n');
fprintf('---------------------------------------\n');
BER_isi = project_4a(OVERSAMPLE, M, k, SNRdB, N, Eb, SNR_linear, nOfIterations, T, taps);

fprintf('\n------------------------------------------------------\n');
fprintf('Discrete LTI Channel with ISI equaliser using FIR filter\n');
fprintf('--------------------------------------------------------\n');
BER_fir = project_4b(OVERSAMPLE, M, k, SNRdB, N, Eb, SNR_linear, nOfIterations, T, taps);

% Plot BER %
textprogressbar('Plot BER: ');
figure(1)
semilogy(SNRdB,BER_isi,'blad:','linewidth',2.5),grid on,hold on;
semilogy(SNRdB,BER_fir,'rd:','linewidth',2.5),grid on,hold on;
title('Discrete LTI BER');
xlabel('SNR(dB)');
ylabel('Bit Error Rate(BER)');
legend('BER ISI', 'BER FIR');
textprogressbar(100, 100);
textprogressbar('done');