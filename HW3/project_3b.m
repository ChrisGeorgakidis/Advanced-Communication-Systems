%*************** Advanced Communication Systems*****************%
%                       CE542, Fall 2018                        %
%                       ECE, UTH, Greece                        %
% File: project_3b.m                                             %
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
rx_antennas = 2;
tx_antennas = 2;
channels = tx_antennas * rx_antennas;

fprintf('---------------------------------------------------\n');
fprintf('| Spatial Multiplexing BER and Goodput Comparison |\n');
fprintf('| Author: Georgakidis Christos                    |\n');
fprintf('| Date: 14/12/2018                                |\n');
fprintf('---------------------------------------------------\n\n');

% All Posible Symbols Generator - kind of look-up table, (look-up vector) %
% each symbol i is in form s(i) = Acos(f(i)) + j * A * sin(f(i)) which is %
% equal to A * exp(j * f(i)), by Euler's formula.                         %
fi = 2.*pi./M.*(0:M-1)';
si = exp(j.*fi);    % all the M possible symbols (complex) %

% TRANSMITTER %
status = 0; % just to show the status %
textprogressbar('Generating Input: ');
for ant = 1:tx_antennas
    % Progress %
    status = status + 1;
    textprogressbar(status, tx_antennas);
    pause(0.1);
    
    % Generate 0,1 with equal probability %
    bits(ant,:) = randi([0 1], 1, N);

    % Create a block of log2(M) bits that will be 'encoded into a symbol %
    splitted_bits = num2str(zeros(ceil(N/log2(M)), log2(M)));

    % Generate the complex signal that will be transmitted %
    for i = 1:log2(M):N
        splitted_bits(ceil(i/log2(M)),:) = num2str(bits(ant, i:i+log2(M)-1));
        symbols(ant, ceil(i/log2(M))) = si(bin2dec(splitted_bits(ceil(i/log2(M)),:)) + 1);
    end
    
    % Oversampling %
    for i = 1:ceil(N/log2(M))
        for o = 1:OVERSAMPLE
            x(ant, i) = symbols(ant, i);
        end
    end
    
    % Calculate the Energy (E) of the symbols %
    E(ant) = sum(abs(x(ant,:)).^2)/(length(x(ant,:)));
end
textprogressbar('done');

% Produce h and noise for each T %
h = zeros(channels,length(x));
w = zeros(tx_antennas,length(x));
for ch = 1:channels
    for packets = 1:(length(x)/T)
        fade = 1 + 1 * (randn(1) + 1i * randn);
        for packet_symbol = 1:T
            h(ch, (packets-1)*T + packet_symbol) = fade;
        end
    end
end
   
BER_mimo = zeros(tx_antennas, length(SNRdB));
H = zeros(rx_antennas, tx_antennas);
y_mimo = zeros(tx_antennas, length(x));
y = zeros(tx_antennas, length(x));

status = 0; % just to show the status %
textprogressbar('Calculating Output: ');
for snr_id = 1:length(SNRdB)
    for mc = 1:nOfIterations    % Monte-Carlo Simulation %
        % Progress %
        status = status + 1;
        textprogressbar(status, length(SNRdB) * nOfIterations);
        pause(0.1);
        
        % CHANNEL %   
        for symbol_id = 1:N
            for rx = 1:rx_antennas
                for tx = 1:tx_antennas
                    H(rx, tx) = h(tx_antennas*(rx - 1) + tx, symbol_id);
                end
            end
            y(:,symbol_id) = awgn(H * x(:,symbol_id), SNRdB(snr_id));
            y_mimo(:,symbol_id) = (H' * H)\H' * y(:,symbol_id);
        end
        
        % RECEIVER %
        for tx = 1:tx_antennas
            samples(tx,:) = y_mimo(tx,OVERSAMPLE:OVERSAMPLE:end);
        end
     
        % DECISION %

        % Find the desired index of the look-up table %
        index = zeros(tx_antennas, ceil(N/log2(M)));
        for tx = 1:tx_antennas
            for i = 1:ceil(N/log2(M))
                [P, I] = min(abs(samples(tx,i) - si));
                index(tx, i) = I;
            end
        end

        % Generate again the block of log2(M) bits %
        % Merge again the blocks of bits in order to construct again the %
        % whole output stream of bits.                                   %
        output = zeros(tx_antennas, ceil(N/log2(M)));
        for tx = 1:tx_antennas
            temp_indx = index(tx,:);
            bits_blocks = dec2bin(temp_indx - 1);
            for i = 1:ceil(N/log2(M))
                for j = 1:log2(M)
                    if bits_blocks(i,j) == '0'
                        
                        output(tx, log2(M)*(i - 1) + j) = 0;
                    else
                        output(tx, log2(M)*(i - 1) + j) = 1;
                    end
                end
            end
        end
       
        % Error Checking %
        error_num = zeros(tx_antennas,1);
        for i = 1:N
            for tx = 1:tx_antennas
                if output(tx, i) ~= bits(tx,i)
                    error_num(tx) = error_num(tx) + 1;
                end
            end
        end

        for tx = 1:tx_antennas
            BER_mimo(tx, snr_id) = BER_mimo(tx, snr_id) + error_num(tx)/N; % MIMO Bit Error Rate %
        end
    end
    
    BER_mimo(:, snr_id) = BER_mimo(:, snr_id)./nOfIterations;
end
textprogressbar('done');

% BER - SNR Graph %
textprogressbar('Plot BER: ');
figure(1)
for tx = 1:tx_antennas
    semilogy(SNRdB,BER_mimo(tx,:),'o','linewidth',2.5),grid on,hold on;
end
title('Spatial Multiplexing BER');
xlabel('SNR(dB)');
ylabel('Bit Error Rate(BER)');
legend('Tx1', 'Tx2');
axis([0 10 10^-5 1]);
textprogressbar(100, 100);
textprogressbar('done');
