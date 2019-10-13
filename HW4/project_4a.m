%*************** Advanced Communication Systems*****************%
%                       CE542, Fall 2018                        %
%                       ECE, UTH, Greece                        %
% File: project_4a.m                                            %
% Authors: Christos Georgakidis (1964)                          %
%***************************************************************%

% clear data from previous run %
% clear
% close all
% clc
% 
% % System Specifications %
% OVERSAMPLE    = 1;
% M             = 2;              % Modulation type %
% k             = log2(M);
% SNRdB         = -100:5:100;           % SNR in dB %
% N             = k*10000;        % number of bits, it is multiplied by log2(M) in%
%                                 % order to ensure that it will be a multiple of it %
% Eb            = 1;              % Energy / bit %
% SNR_linear    = 10.^(SNRdB/10); % convert SNR from db to linear: SNR(dB) = 10log10(SNR(linear)) %
% nOfIterations = 10;             % num of iterations for Monte Carlo Simulation %
% T = 5 * k; % The time the channels remains the same %
% taps = 5;   % Number of taps %
% 
% 
% fprintf('----------------------------------\n');
% fprintf('| Discrete LTI Channel with ISI  |\n');
% fprintf('| Author: Georgakidis Christos   |\n');
% fprintf('| Date: 20/12/2018               |\n');
% fprintf('----------------------------------\n\n');
% 

function BER = project_4a(OVERSAMPLE, M, k, SNRdB, N, Eb, SNR_linear, nOfIterations, T, taps)

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
x = zeros(ceil(N/log2(M)) * OVERSAMPLE, 1);
for i = 1:ceil(N/log2(M))
    for o = 1:OVERSAMPLE
        x((i-1) + o) = symbols(i);
    end
end

% Calculate the Energy (E) of the symbols %
E = sum(abs(x).^2)/(length(x));

textprogressbar(100, 100);
pause(0.1);
textprogressbar('done');

textprogressbar('Generating Channel: ');

% Construct the h vectors %
H = 1 + 1 * (randn(1,ceil((N + taps - 1)/T) * T) + 1i*randn(1,ceil((N + taps - 1)/T) * T));

textprogressbar(100, 100);
pause(0.1);
textprogressbar('done');

status = 0;
textprogressbar('Calculating Output: ');
BER = zeros(1, length(SNRdB));
for snr_id = 1:length(SNRdB)
    for mc = 1:nOfIterations    %Monte-Carlo Simulation %
        % Progress %
        status = status + 1;
        textprogressbar(status, length(SNRdB) * nOfIterations);
        pause(0.1);

        % CHANNEL %
        y_isi = channel_effect(H, x, length(x), taps, T);
        
        y_temp = awgn(y_isi, SNRdB(snr_id));
        
        % RECEIVER %
        
        % Equalisation %
        y = equaliser(H, y_temp, length(x), taps, T);
        
        % Sampling %
        samples = y(OVERSAMPLE:OVERSAMPLE:end);
        
        % DECISION %
        
        % Find the desired index of the look-up table %
        index = [];
        for idx = 1:ceil(N/log2(M))
            [P, I] = min(abs(samples(idx) - si));
            index = [index I];
        end
        
        % Generate again the block of log2(M) bits %
        bits_blocks = dec2bin(index - 1);
        
        % Merge again the blocks of bits in order to construct again the %
        % whole output stream of bits.                                   %
        output = [];
        for i = 1:ceil(N/log2(M))
            for jj = 1:log2(M)                
                if bits_blocks(i, jj) == '0'
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
        
        BER(snr_id) = BER(snr_id) + error_num/N;
       
    end
    
    BER(snr_id) = BER(snr_id)/nOfIterations;
        
end
textprogressbar('done');

% % BER - SNR Graph %
% textprogressbar('Plot BER: ');
% 
% % Rayleigh Theoretical BER
% E=sum(abs(x).^2)/(length(x));
% No=E./SNR_linear;
% G = 2*E./No;
% BER_Rayleigh_theor = (1/2)*(1-sqrt(G./(G+1)));
% 
% figure(2)
% semilogy(SNRdB,BER,'blad:','linewidth',2.5),grid on,hold on;
% % semilogy(SNRdB,BER_Rayleigh_theor, 'blad-','linewidth',2.5); 
% title('Discrete LTI BER');
% xlabel('SNR(dB)');
% ylabel('Bit Error Rate(BER)');
% % legend('MRC', 'Rayleigh-theoretical');
% % axis([0 length(SNRdB) 10^-5 1]);
% textprogressbar(100, 100);
% textprogressbar('done');
end

%%% Channel Effect %%%
% This function creates the channel's effect to the symbols applying ISI %
function y = channel_effect(H, x, N, taps, T)
    y = zeros(1,length(x) + taps - 1);
    
    % extend x with zeros s.t: 
    %     x_exp = [ 0 0 0 0 x(1) x(2) ... x(N) 0 0 0 0 ]
    % we do that in order to use the same formula with all symbol bits
    temp = zeros(taps - 1, 1);
    x_exp = [ temp; x; temp ];
    
    offset = 1;
    for k = 1:(N + taps - 1)
        if (mod(k - taps - 1, T) == 0) && (k > T)
            offset = offset + taps;
        end
        y(k) = H(offset:offset + taps - 1) * (x_exp(k + taps - 1:-1:k));
    end
end

%%% Equaliser %%%
% This function performs equalisation using the least-square method %
function y_eq = equaliser(H, y, N, taps, T)
    y_eq = zeros(1,N);
    
    offset = 1;
    for k = 1:N
        if (mod(k - taps - 1, T)==0) && (k > T)
            offset = offset + taps;
        end
        for s = 0:4
            y_eq(k) = y_eq(k) + H(offset+s)' * y(k+s) / norm(H(offset+s));
        end
        y_eq(k) = mean(y_eq(k));
    end
end
