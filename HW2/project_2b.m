%*************** Advanced Communication Systems*****************%
%                       CE542, Fall 2018                        %
%                       ECE, UTH, Greece                        %
% File: project_2b.m                                             %
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
SNRdB         = 0:20;           % SNR in dB %
N             = k*10000;        % number of bits, it is multiplied by log2(M) in%
                                % order to ensure that it will be a multiple of it %
Eb            = 1;              % Energy / bit %
SNR_linear    = 10.^(SNRdB/10); % convert SNR from db to linear: SNR(dB) = 10log10(SNR(linear)) %
nOfIterations = 10;             % num of iterations for Monte Carlo Simulation %
T = 1 * k;

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
E = sum(abs(x).^2)/(length(x)); %Calculate actual symbol energy

% Produce h and noise for each T %
h = zeros(1,length(x));
w = zeros(1,length(x));
for i = 1:(length(x)/T)
    fadeI = 1 + sqrt(1)*randn(1);
    fadeQ = 1 + sqrt(1)*randn(1);
    fade = 1 + 1 * ( randn(1) + 1i * randn);
    for ii = 1:T
        h((i-1)*T + ii) = fade;
    end
end
    
BER = zeros(1, length(SNRdB));
BER_fade = zeros(1, length(SNRdB));

for c = 0:10
    BER_fadeInv(c+1,:) = zeros(1, length(SNRdB));
    for ii = 1:length(SNRdB)
    
        No=E/SNR_linear(ii); %Find the noise spectral density
        w = sqrt(No/2)*(randn(1,length(x))+1i*randn(1,length(x)));%computed noise  

        for mc = 1:nOfIterations    % Monte-Carlo Simulation %
            %  Channel  %
            y = awgn(x, SNRdB(ii));

            y_fade = h.*x + w; % No Channel Inversion %
            % y_fadeInv = y_fade./h; % Channel Inversion %
            ci = (h+c)'/norm(h+c);
            for l = 1:length(y_fade)
                y_fadeInv(l) = y_fade(l) * ci(l);
            end

            % I commented out the construction of constellation diagram because %
            % for many input bits and monte carlo iterations, my laptop was     %
            % killing matlab because it was frozen. In order to check them out  %
            % remove it from comment.                                           %

            %scatterplot(y);

            % see how channel affect the signal %
            if (SNRdB(ii) == 4 && mc == 1 && c == 1)
                y_plot = zeros(1, T);
                y_plot_inv = zeros(1, T);
                x_plot = zeros(1, T);

                for i = 1:T
                    y_plot(i) = y_fade(i);
                    x_plot(i) = x(i);
                end

                for i = 1:T
                    y_plot_inv(i) = y_fadeInv(i);
                end

                figure(1);
                scatter(real(x_plot), real(y_plot), 'r');
                hold on;
                scatter(real(x_plot), real(y_plot_inv), 'b');
                grid on;

            end

            %  Receiver  %
            samples = y(OVERSAMPLE:OVERSAMPLE:end); % sampling %

            fade_samples = y_fade(OVERSAMPLE:OVERSAMPLE:end); % Rayleigh samples %
            fadeInv_samples = y_fadeInv(OVERSAMPLE:OVERSAMPLE:end); % Rayleigh samples with channel inversion %

            % Decision %

            % Find the desired index of the look-up table %
            index = [];
            index_fade = [];
            index_fadeInv = [];
            for i = 1:ceil(N/log2(M))
                [P, I] = min(abs(samples(i) - si));
                [Pfade, Ifade] = min(abs(fade_samples(i) - si));
                [Pfade_inv, Ifade_inv] = min(abs(fadeInv_samples(i) - si));
                index = [index I];
                index_fade = [index_fade Ifade];
                index_fadeInv = [index_fadeInv Ifade_inv];
            end

            % Generate again the block of log2(M) bits %
            bits_blocks = dec2bin(index - 1);
            bits_blocks_fade = dec2bin(index_fade - 1);
            bits_blocks_fadeInv = dec2bin(index_fadeInv - 1);

            % Merge again the blocks of bits in order to construct again the %
            % whole output stream of bits.                                   %
            output = [];
            output_fade = [];
            output_fadeInv = [];
            for i = 1:ceil(N/log2(M))
                for j = 1:log2(M)
                    if bits_blocks(i, j) == '0'
                        output = [output 0];
                    else
                        output = [output 1];
                    end

                    if bits_blocks_fade(i, j) == '0'
                        output_fade = [output_fade 0];
                    else
                        output_fade = [output_fade 1];
                    end

                    if bits_blocks_fadeInv(i, j) == '0'
                        output_fadeInv = [output_fadeInv 0];
                    else
                        output_fadeInv = [output_fadeInv 1];
                    end
                end
            end

            % Error Checking %
            error_num = 0;
            error_num_fade = 0;
            error_num_fadeInv = 0;
            for i = 1:N
                if output(i) ~= bits(i)
                    error_num = error_num + 1;
                end

                if output_fade(i) ~= bits(i)
                    error_num_fade = error_num_fade + 1;
                end

                if output_fadeInv(i) ~= bits(i)
                    error_num_fadeInv = error_num_fadeInv + 1;
                end
            end

            BER(ii) = BER(ii) + error_num/N; % Bit Error Rate %
            BER_fade(ii) = BER_fade(ii) + error_num_fade/N;
            BER_fadeInv(c+1,ii) = BER_fadeInv(c+1,ii) + error_num_fadeInv/N;
        end
        BER(ii) = BER(ii)/nOfIterations;
        BER_fade(ii) = BER_fade(ii)/nOfIterations;
        BER_fadeInv(c+1,ii) = BER_fadeInv(c+1,ii)/nOfIterations;
    end
end

% Rayleigh Theoretical BER
E=sum(abs(x).^2)/(length(x));
No=E./SNR_linear;
G = 2*E./No;
BER_Rayleigh_theor = (1/2)*(1-sqrt(G./(G+1)));

% AWGN Theoretical BER %
BER_theor = (1/2)*erfc(sqrt(SNR_linear));

% BER - SNR Graph %
figure(2)
semilogy(SNRdB,BER,'o','linewidth',2.5),grid on,hold on;
semilogy(SNRdB,BER_fade,'o','linewidth',2.5),grid on,hold on;
for c = 0:10
    semilogy(SNRdB,BER_fadeInv(c+1,:),'o','linewidth',2.5),grid on,hold on;
end
semilogy(SNRdB,BER_theor,'-*r','linewidth',2.5);
semilogy(SNRdB,BER_Rayleigh_theor, 'blad-','linewidth',2.5); 
title('Bit Error Rate vs SNR for mPSK');
xlabel('SNR(dB)');
ylabel('Bit Error Rate(BER)');
legend('AWGN', 'Rayleigh', 'Channel Inversion', 'AWGN-theoretical', 'Rayleigh-theoretical');
axis([0 10 10^-5 1]);

% Note %
% For some reason, I don't know the simulation is below the theoretical %
% Rayleigh line.                                                        %

% Task B2 (I Didn't have time to implement it in code) %
% One way to improve the BER is the following: Every time the channel    %
% changes the transmitter can send to receiver some predefined symbols.  %
% Doing that, the receiver knowing the expected value of y can solve     %
% the equation y = hmod * x + w and find an average value of h, since we %
% have the unknown w.                                                    %


