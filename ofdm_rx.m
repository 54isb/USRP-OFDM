function data_est = ofdm_rx(rx, Nc, Cp, N)
% rx : received signal 
% Nc : Num. of subcarriers (must-be-even)
% M  : Num. of OFDM Symbols
% data : Nc × M QAM Symbol Matrix in Freq. Domain
% Cp : Cyclic prefix length
% N  : Num. of IFFT points

    %Precheck
    assert(mod(Nc,2)==0, 'Nc must be even.');
    assert(size(rx,1)==N+Cp, 'Col. of input matrix must be N+Cp.');
    assert(mod(N-Nc,2)==0, '(N-Nc) must be even.');

    %Obtain Num. of OFDM symbols
    M = size(rx,2);

    %Delete CP -> N x M 
    y = rx(Cp+1:end, :, :);        

    %FFT: 1/sqrt(N/Nc) is normalization for oversampling
    Y = sqrt(Nc) * 1/N * fft(y, N, 1); 

    %Get corresponding subcarriers 
    half = Nc/2;
    data_est = [Y(1:half,:,:);Y(end-half+1:end,:,:)];

end