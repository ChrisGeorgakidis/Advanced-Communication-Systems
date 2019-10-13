%*************** Advanced Communication Systems*****************%
%                       CE542, Fall 2018                        %
%                       ECE, UTH, Greece                        %
% File: project_3c.m                                             %
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

fprintf('----------------------------------------------------------\n');
fprintf('| MRC vs Spatial Multiplexing BER and Goodput Comparison |\n');
fprintf('| Author: Georgakidis Christos                           |\n');
fprintf('| Date: 14/12/2018                                       |\n');
fprintf('----------------------------------------------------------\n\n');

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
    
    % Generate 0, 1 with equal probability %
    bits(ant,:) = randi([0 1], 1, N);

    % Create a block of log2(M) bits that will be 'encoded' into a symbol %
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

h_mrc = zeros(rx_antennas, length(x));
offset = 0;
for rx = 1:rx_antennas
    h_mrc(rx,:) = h(offset + rx,:);
    offset = offset + rx_antennas;
end

BER_mimo = zeros(tx_antennas, length(SNRdB));
BER_mrc = zeros(1, length(SNRdB));
goodput_sm = zeros(tx_antennas, length(SNRdB));
goodput_mrc = zeros(1, length(SNRdB));
y_mimo = zeros(tx_antennas, length(x));
y_mrc = zeros(1, length(x));
y1 = zeros(tx_antennas, length(x));
y2 = zeros(1, length(x));
H = zeros(rx_antennas, tx_antennas);
h_sum = zeros(1, length(x));
h_norm = zeros(1, length(x));

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

            % Spatial Multiplexing %
            for rx = 1:rx_antennas
                for tx = 1:tx_antennas
                    H(rx, tx) = h(tx_antennas*(rx - 1) + tx, symbol_id);
                end
            end

            y1(:,symbol_id) = awgn(H * x(:,symbol_id), SNRdB(snr_id));
            y_mimo(:,symbol_id) = (H' * H)\H' * y1(:,symbol_id);

            % MRC %
            for rx = 1:rx_antennas
                h_sum(symbol_id) = h_sum(symbol_id) + norm(h_mrc(rx, symbol_id));
            end
            h_norm(symbol_id) = norm(h_mrc(:, symbol_id));
        end
        y2 = awgn(x(1,:), SNRdB(snr_id));
        y_mrc = h_sum/h_norm * y2;

        % RECEIVER %

        % spatial multiplexing %
        for tx = 1:tx_antennas
            mimo_samples(tx,:) = y_mimo(tx,OVERSAMPLE:OVERSAMPLE:end);
        end

        % mrc %
        mrc_samples = y_mrc(OVERSAMPLE:OVERSAMPLE:end);


        % DECISION %

        % Find the desired index of the look-up table %

        % spatial multiplexing %
        index_mimo = zeros(tx_antennas, ceil(N/log2(M)));
        for tx = 1:tx_antennas
            for i = 1:ceil(N/log2(M))
                [P_mimo, I_mimo] = min(abs(mimo_samples(tx,i) - si));
                index_mimo(tx, i) = I_mimo;
            end
        end

        % Generate again the block of log2(M) bits %
        % Merge again the blocks of bits in order to construct again the %
        % whole output stream of bits.                                   %
        output_mimo = zeros(tx_antennas, ceil(N/log2(M)));
        for tx = 1:tx_antennas
            temp_indx = index_mimo(tx,:);
            bits_blocks_mimo = dec2bin(temp_indx - 1);
            for i = 1:ceil(N/log2(M))
                for j = 1:log2(M)
                    if bits_blocks_mimo(i,j) == '0'

                        output_mimo(tx, log2(M)*(i - 1) + j) = 0;
                    else
                        output_mimo(tx, log2(M)*(i - 1) + j) = 1;
                    end
                end
            end
        end

        % Error Checking %
        error_num_mimo = zeros(tx_antennas,1);
        for tx = 1:tx_antennas
            for packets = 1:(length(x)/T)
                error_packet_mimo = 0;
                for packet_symbol = 1:T
                    if output_mimo(tx, (packets-1)*T + packet_symbol) ~= bits(tx,(packets-1)*T + packet_symbol)
                        error_num_mimo(tx) = error_num_mimo(tx) + 1;
                        error_packet_mimo = 1;
                    end
                end

                if error_packet_mimo ~= 1
                    goodput_sm(tx, snr_id) = goodput_sm(tx, snr_id) + 1;
                end
            end
        end

        for tx = 1:tx_antennas
            BER_mimo(tx, snr_id) = BER_mimo(tx, snr_id) + error_num_mimo(tx)/N; % MIMO Bit Error Rate %
        end

        % mrc %
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
        for packets = 1: (length(x)/T)
            error_packet_mrc = 0;
            for packet_symbol = 1:T
                if output_mrc((packets-1)*T + packet_symbol) ~= bits(1,(packets-1)*T + packet_symbol)
                    error_num_mrc = error_num_mrc + 1;
                    error_packet_mrc = 1;
                end
            end

            if error_packet_mrc ~= 1
                goodput_mrc(snr_id) = goodput_mrc(snr_id) + 1;
            end
        end

        BER_mrc(snr_id) = BER_mrc(snr_id) + error_num_mrc/N;
    end

    BER_mimo(:, snr_id) = BER_mimo(:, snr_id)./nOfIterations;
    goodput_sm(:, snr_id) = goodput_sm(:, snr_id)./nOfIterations;

    BER_mrc(snr_id) = BER_mrc(snr_id)/nOfIterations;
    goodput_mrc(snr_id) = goodput_mrc(snr_id)/nOfIterations;

end
textprogressbar('done');

% Goodput - SNR Graph %
textprogressbar('Plot Goodput: ');
figure(1)
for tx = 1:tx_antennas
    semilogy(SNRdB,goodput_sm(tx,:),'o','linewidth',2.5),grid on,hold on;
end
semilogy(SNRdB,goodput_mrc,'*','linewidth',2.5),grid on,hold on;
title('MRC vs Spatial Multiplexing Goodput');
xlabel('SNR(dB)');
ylabel('Goodput');
legend('SM-Tx1','SM-Tx2', 'MRC');
textprogressbar(100, 100);
textprogressbar('done');


% BER - SNR Graph %
textprogressbar('Plot BER: ');
figure(2)
for tx = 1:tx_antennas
    semilogy(SNRdB,BER_mimo(tx,:),'o','linewidth',2.5),grid on,hold on;
end
semilogy(SNRdB,BER_mrc,'*','linewidth',2.5),grid on,hold on;
title('MRC vs Spatial Multiplexing BER');
xlabel('SNR(dB)');
ylabel('Bit Error Rate(BER)');
legend('SM-Tx1','SM-Tx2', 'MRC');
axis([0 10 10^-5 1]);
textprogressbar(100, 100);
textprogressbar('done');

