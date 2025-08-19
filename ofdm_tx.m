function tx = ofdm_tx(data, Nc, Cp, N)
% OFDM Transmitter
% Nc : Num. of subcarriers (must-be-even)
% M  : Num. of OFDM Symbols
% data : Nc × M QAM Symbol Matrix in Freq. Domain
% Cp : Cyclic prefix length
% N  : Num. of IFFT points
% Output becomes the sequence with the length (N+Cp)xM in time domain   
    assert(mod(Nc,2)==0, 'Nc must be even.');
    [nRow, M, T] = size(data);
    assert(nRow==Nc, 'data の行数 must be Nc.');
    assert(mod(N-Nc,2)==0, '(N-Nc) must be even.');

    %subcarrier allocation
    half = Nc/2;
    X = [data(1:half,:,:) ; zeros(N-Nc, M, T) ; data(half+1:end,:,:)];

    %IFFT Operation (note: sqrt(N/Nc) is normalization for oversampling)
    x = N * 1/sqrt(Nc) * ifft(X, N); 

    %Add cyclic prefix
    cp_block = x(end-Cp+1:end, :,:);  
    tx = [cp_block; x];              

end