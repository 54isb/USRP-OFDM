clear;clc;close all;

%OFDM parameter setting ******************************************
Nc = 140;    % Num. of subcarriers
Cp = 70;    % CP length
N  = 512;   % IFFT size
M  = 5;    % Num. of OFDM symbols
%PHY frame setting **********************************************
Mp = 4;           % Modulation order (M-ary QAM)
k = log2(Mp);     % Bits per symbol
numBits = k*Nc; % Number of bits to process
preamlen = 139;   % Length of ZF sequence
datalen = numBits/k; % Length of payload
ZCroot = 25; %Root of ZC sequence
%USRP Setting **********************************************
use_usrp = 0; %0: USRP OFF (Debug) 1: USRP ON
Fc = 5.0e9;   %Center frequency
TxGain = 60;     % Transmission Gain (0 dB to 60 dB for N320)
RxGain = 60;     % Receiver Gain (0 dB to 60 dB for N320)
sps = 125e+6;    % Default sampling speed of N320 (196,851 Hz to 125 MHz, 200.00 MHz, 245.76 MHz, 250.00 MHz (default))

%Output Setting
fprintf('Bandwidth: %.4f [MHz], Frame Duration %.4f [usec] \n', (sps/N * Nc)/1e+6, (1/sps * (N+Nc)*(M+2))/1e-6);

if use_usrp == 1
    %%USRP setting ************************************************
    bbtrx = basebandTransceiver("My USRP N320",CaptureDataType="double",TransmitDataType="double",TransmitRadioGain=TxGain,CaptureRadioGain=RxGain,CaptureAntennas="RF0:RX2",TransmitAntennas="RF0:TX/RX",TransmitCenterFrequency=Fc,CaptureCenterFrequency=Fc);
    %%USRP setting ************************************************
    [pre_signal,~] = capture(bbtrx,milliseconds(2)); %Received signal stream from USRP
    %Obtain inband noise
    noise_rec = reshape(pre_signal(1:(N+Cp)*10*M),N+Cp,10*M); %Reshape for FFT
    inband_noise = ofdm_rx(noise_rec, Nc, Cp, N); %Get inband noise at the subcarriers
    noise_var = rms(inband_noise(:))^2; %Calculate noise power
end

% Training and Pilot Sequences
preamble = [zadoffChuSeq(ZCroot,preamlen);0]; %ZC sequence 1 for frame sync.
preamble2 = [zadoffChuSeq(ZCroot,preamlen-2);0;0;0]; %ZC sequence 2 for frame sync.
preamble3 = [zadoffChuSeq(ZCroot,preamlen-6);0;0;0;0;0;0;0]; %ZC sequence 3 for frame sync. 
pilot = 2 * randi([0, 1],Nc,1) - 1; %all one pilot symbols with dithering for PAPR reduction

% Generate QAM symbols
dataIn = randi([0 1],numBits,M); % Generate vector of binary data
dataSymbolsIn = bit2int(dataIn,k); % Binary -> GF(Mp)
payload = qammod(dataSymbolsIn, Mp, UnitAveragePower=true); %data modulation

%make PHY frame (Nc x (sync 3 + pilot 1 + data M)
frame = [preamble,preamble2,preamble3,pilot,payload]; 

% OFDM Symbol Modulation
tx_signal_matrix = ofdm_tx(frame, Nc, Cp, N); %ofdm modulation
tx_signal = tx_signal_matrix(:); %parallel to serial
%input waveform (to USRP) must be in range [-1, 1]
tx_signal(:) = 1/max(abs(tx_signal(:))) * tx_signal(:);
%tx_signal(1:(N+Cp)*3) = 1/abs(max(tx_signal(1:(N+Cp)*3))) * tx_signal(1:(N+Cp)*3); % Sync part 
%tx_signal((N+Cp)*3+1:end) = 1/abs(max(tx_signal(1:(N+Cp)*3))) * tx_signal((N+Cp)*3+1:end); % Pilot + Data part

%USPR Operation Begin*****************************************************
if use_usrp == 1
    transmit(bbtrx,tx_signal,"continuous"); %Transmit from USRP
    [rx_signal,~] = capture(bbtrx,milliseconds(2)); %Received signal stream from USRP
    stopTransmission(bbtrx); %Close connection of USRP
else
    %Reception for debug    
    rx_signal = tx_signal;
end
%USPR Operation End*******************************************************

%Frame synchronization based on ZC sequence
syncsize = ceil(size(rx_signal,1)/2); %Search up to half the length of the received frame
[val, pos] = max(abs(xcorr(rx_signal,reshape(tx_signal_matrix(:,1:3),[],1),syncsize))); %Preamble search
pos = pos - syncsize; %Eliminate the effect of the negative lags of xcorr

%Obtain the payload part (w/ M OFDM symbols)
sync_rx_signal=rx_signal(pos+(N+Cp)*3:pos+(N+Cp)*(M+4)-1);

%Reshape (N+Cp)(M+1) -> (N+Cp) x (M+1)
sync_rx_signal_matrix = reshape(sync_rx_signal, N+Cp, M+1);

%OFDM Symbol Demodulation
rec_symbols = ofdm_rx(sync_rx_signal_matrix, Nc, Cp, N);

%ZF Equlization using all-one pilots
equalized_symbols = rec_symbols(:,2:end) ./ (rec_symbols(:,1) .* pilot); 
demod_data = qamdemod(equalized_symbols,Mp,UnitAveragePower=true);

%decimal -> binary
dataOut = int2bit(demod_data,k);
%Calculate the number of errors
[numErrors,ber] = biterr(dataIn,dataOut);
%Output BER and the number of errors
if use_usrp == 1
    scatterplot(equalized_symbols(:));
    figure(2);periodogram(rx_signal,'centered');
end
fprintf('The bit error rate is %5.2e (%d errors).\n',ber,numErrors)