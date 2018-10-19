%*************** Advanced Communication Systems*****************%
%                       CE542, Fall 2018                        %
%                       ECE, UTH, Greece                        %
% File: project_1.m                                             %
% Authors: Christos Georgakidis (1964)                          %
%***************************************************************%

% clear data from previous run %
clear
close all
clc

% System Specifications %
OVERSAMPLE    = 4;
M             = 2;              % Modulation type %
k             = log2(M);
SNRdB         = 0:20;           % SNR in dB %
N             = k*10000;        % number of bits, it is multiplied by log2(M) in%
                                % order to ensure that it will be a multiple of it %
Eb            = 1;              % Energy / bit %
SNR_linear    = 10.^(SNRdB/10); % convert SNR from db to linear: SNR(dB) = 10log10(SNR(linear)) %
nOfIterations = 10;             % num of iterations for Monte Carlo Simulation %


% All Posible Symbols Generator - kind of look-up table, (look-up vector) %
% each symbol i is in form s(i) = Acos(f(i)) + j * A * sin(f(i)) which is %
% equal to A * exp(j * f(i)), by Euler's formula.                         %
fi = 2.*pi./M.*(0:M-1)';
si = exp(j.*fi);    % all the M possible symbols (complex) %

%  Transmitter  %
bits = randi([0 1], 1, N);  % generating 0, 1 with equal probability %

splitted_bits = num2str(zeros(ceil(N/log2(M)), log2(M)));   % create a block of log2(M) bits that will be %
                                                            % "encoded" into a symbol                     %

% Generates the complex signal that will be transmitted %
for i = 1:log2(M):N
    splitted_bits(ceil(i/log2(M)),:) = num2str(bits(i:i+log2(M)-1));
    symbols(ceil(i/log2(M))) = si(bin2dec(splitted_bits(ceil(i/log2(M)),:)) + 1);
end

% The complex signal vector is goind to be oversampled %
x = [];
for i = 1:ceil(N/log2(M))
    for o = 1:OVERSAMPLE
        x = [x symbols(i)];
    end
end

BER = zeros(1, length(SNRdB));
for ii = 1:length(SNRdB)
    for mc = 1:nOfIterations    % Monte-Carlo Simulation %
        %  Channel  %
        y = awgn(x, SNRdB(ii));

        % I commented out the construction of constellation diagram because %
        % for many input bits and monte carlo iterations, my laptop was     %
        % killing matlab because it was frozen. In order to check them out  %
        % remove it from comment.                                           %

        %scatterplot(y);

        %  Receiver  %
        samples = y(OVERSAMPLE:OVERSAMPLE:end); % sampling %

        % Decision %

        % Find the desired index of the look-up table %
        index = [];
        for i = 1:ceil(N/log2(M))
            [P, I] = min(abs(samples(i) - si));
            index = [index I];
        end

        % Generate again the block of log2(M) bits %
        bits_blocks = dec2bin(index - 1);

        % Merge again the blocks of bits in order to construct again the %
        % whole output stream of bits.                                   %
        output = [];
        for i = 1:ceil(N/log2(M))
            for j = 1:log2(M)
                if bits_blocks(i, j) == '0'
                    output = [output 0];
                else
                    output = [output 1];
                end
            end
        end

        % Error Checking %
        error_num = 0;
        for i = 1:N
            if output(i) ~= bits(i)
                error_num = error_num + 1;
            end
        end

        BER(ii) = BER(ii) + error_num/N; % Bit Error Rate %
    end
    BER(ii) = BER(ii)/nOfIterations;
end

BER_theor = (1/2)*erfc(sqrt(SNR_linear));

% BER - SNR Graph %
figure(1)
semilogy(SNRdB,BER,'o','linewidth',2.5),grid on,hold on;
semilogy(SNRdB,BER_theor,'r','linewidth',2.5);
title('Bit Error Rate vs SNR for mPSK');
xlabel('SNR(dB)');
ylabel('Bit Error Rate(BER)');
legend('simulation','theorytical');
axis([0 10 10^-5 1]);