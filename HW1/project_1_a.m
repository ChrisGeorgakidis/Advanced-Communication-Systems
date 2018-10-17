%*************** Advanced Communication Systems*****************%
%                       CE542, Fall 2018                        %
%                       ECE, UTH, Greece                        %
% File: project_1                                               %
% Authors: Christos Georgakidis (1964)                          %
%***************************************************************%

% pkg load communications
clear all
close all

SNRdB = 0:20; % SNR in dB %
N = 100000; % number of bits %

% System Specifications %
Eb = 1; % Energy / bit %

SNR_linear = 10.^(SNRdB/10);
nOfIterations = 1;  % num of iterations for Monte Carlo Simulation %

%  Transmitter  %
bits = rand(1, N) > 0.5; % generating 0, 1 with equal probability %
s = 2*bits - 1; % BPSK modulation 0 -> -1, 1 -> 1 %
x = complex(s);

for ii = 1:length(SNRdB)
    %  Channel  %
    y = awgn(x, SNRdB(ii));

%     scatterplot(y);

    %  Receiver  %
    output = real(y)>0; % decision %

    % counting the errors %
    error_num = 0;
    for i = 1:N
        if output(i) ~= bits(i)
            error_num = error_num + 1;
        end
    end
    BER(ii) = error_num/N; % Bit Error Rate %
end

BER_exp = (1/2)*erfc(sqrt(SNR_linear));

% BER - SNR Graph %
figure(1)
%plot start
semilogy(SNRdB,BER,'o','linewidth',2.5),grid on,hold on;
semilogy(SNRdB,BER_exp,'r','linewidth',2.5);
title(' curve for Bit Error Rate verses  SNR for Binary PSK modulation');
xlabel(' SNR(dB)');
ylabel('BER');
legend('simulation','theorytical');
axis([0 10 10^-5 1]);





