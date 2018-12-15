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
OVERSAMPLE    = 1;
M             = 2;              % Modulation type %
k             = log2(M);
SNRdB         = 0:10;           % SNR in dB %
N             = k*10000;        % number of bits, it is multiplied by log2(M) in%
                                % order to ensure that it will be a multiple of it %
Eb            = 1;              % Energy / bit %
SNR_linear    = 10.^(SNRdB/10); % convert SNR from db to linear: SNR(dB) = 10log10(SNR(linear)) %
nOfIterations = 10;             % num of iterations for Monte Carlo Simulation %
T = 10 * k;
diversity_branches = 2;

fprintf('----------------------------------\n');
fprintf('| MRC BER and Goodput Comparison |\n');
fprintf('| Author: Georgakidis Christos   |\n');
fprintf('| Date: 14/12/2018               |\n');
fprintf('----------------------------------\n\n');


% All Posible Symbols Generator - kind of look-up table, (look-up vector) %
% each symbol i is in form s(i) = Acos(f(i)) + j * A * sin(f(i)) which is %
% equal to A * exp(j * f(i)), by Euler's formula.                         %
fi = 2.*pi./M.*(0:M-1)';
si = exp(j.*fi);    % all the M possible symbols (complex) %

%  TRANSMITTER  %
textprogressbar('Generating Input: ');

% Generate 0,1 with equal probability %
bits = randi([0 1], 1, N);

% Create a block of log2(M) bits that will be 'encoded into a symbol %
splitted_bits = num2str(zeros(ceil(N/log2(M)), log2(M)));

% Generate the complex signal that will be transmitted %
for i = 1:log2(M):N
    splitted_bits(ceil(i/log2(M)),:) = num2str(bits(i:i+log2(M)-1));
    symbols(ceil(i/log2(M))) = si(bin2dec(splitted_bits(ceil(i/log2(M)),:)) + 1);
end

% Oversampling %
x = [];
for i = 1:ceil(N/log2(M))
    for o = 1:OVERSAMPLE
        x = [x symbols(i)];
    end
end

% Calculate the Energy (E) of the symbols %
E = sum(abs(x).^2)/(length(x));

textprogressbar(100, 100);
pause(0.1);
textprogressbar('done');

% Produce h and noise for each T %
h = zeros(diversity_branches,length(x));
w = zeros(1,length(x));
for d = 1:diversity_branches
    for i = 1:(length(x)/T)
        fade = 1 + 1 * (randn(1) + 1i * randn);
        for snr_id = 1:T
            h(d, (i-1)*T + snr_id) = fade;
        end
    end
end

BER_mrc = zeros(1, length(SNRdB));
y = zeros(1, length(x));
y_mrc = zeros(1, length(x));
h_sum = zeros(1, length(x));
h_norm = zeros(1, length(x));

status = 0;
textprogressbar('Calculating Output: ');
for snr_id = 1:length(SNRdB)
    
    No=E/SNR_linear(snr_id); %Find the noise spectral density
    w = sqrt(No/2)*(randn(1,length(x))+1i*randn(1,length(x)));%computed noise  

    for mc = 1:nOfIterations    % Monte-Carlo Simulation %
        % Progress %
        status = status + 1;
        textprogressbar(status, length(SNRdB) * nOfIterations);
        pause(0.1);
        
        % CHANNEL % 
        for symbol_id = 1:N
            h_sum(symbol_id) = 0;
            for d = 1:diversity_branches
                h_sum(symbol_id) = h_sum(symbol_id) + norm(h(d,symbol_id));
            end
            h_norm(symbol_id) = norm(h(:, symbol_id));
        end
        y = awgn(x, SNRdB(snr_id));
        y_mrc = h_sum/h_norm * y;
        
        
        % RECEIVER %
        
        mrc_samples = y_mrc(OVERSAMPLE:OVERSAMPLE:end); % MRC samples % 
        
        % DECISION %

        % Find the desired index of the look-up table %
        index_mrc = [];
        for i = 1:ceil(N/log2(M))
            [Pmrc, Imrc] = min(abs(mrc_samples(i) - si));
            index_mrc = [index_mrc Imrc];
        end

        % Generate again the block of log2(M) bits %
        bits_blocks_mrc = dec2bin(index_mrc - 1);
       
        % Merge again the blocks of bits in order to construct again the %
        % whole output stream of bits.                                   %
        output_mrc = [];
        for i = 1:ceil(N/log2(M))
            for j = 1:log2(M)                
                if bits_blocks_mrc(i, j) == '0'
                    output_mrc = [output_mrc 0];
                else
                    output_mrc = [output_mrc 1];
                end
            end
        end

        % Error Checking %
        error_num_mrc = 0;
        for i = 1:N            
            if output_mrc(i) ~= bits(i)
                error_num_mrc = error_num_mrc + 1;
            end
        end
        
        BER_mrc(snr_id) = BER_mrc(snr_id) + error_num_mrc/N;
    end

    BER_mrc(snr_id) = BER_mrc(snr_id)/nOfIterations;
end
textprogressbar('done');

% BER - SNR Graph %
textprogressbar('Plot BER: ');

% Rayleigh Theoretical BER
E=sum(abs(x).^2)/(length(x));
No=E./SNR_linear;
G = 2*E./No;
BER_Rayleigh_theor = (1/2)*(1-sqrt(G./(G+1)));

figure(2)
semilogy(SNRdB,BER_mrc,'*','linewidth',2.5),grid on,hold on;
semilogy(SNRdB,BER_Rayleigh_theor, 'blad-','linewidth',2.5); 
title('MRC BER');
xlabel('SNR(dB)');
ylabel('Bit Error Rate(BER)');
legend('MRC', 'Rayleigh-theoretical');
axis([0 10 10^-5 1]);
textprogressbar(100, 100);
textprogressbar('done');
