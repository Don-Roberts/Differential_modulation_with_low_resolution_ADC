function [decodedSym,shat_mag]=ML1BitAlamouti(Y,NSub, Nt,H_est,ModulationType, SNR)
    %-1 = 1;
    N = size(Y,1);
    U = size(Y,2);

    NSub = N/2;
    shat = zeros(1,2);
    shat_mag = zeros(NSub,1);
    m= 1;


    YSignreal = sign(real(Y));
    YSignimg = sign(imag(Y));
    
    %YSignreal = (real(Y));
    %YSignimg = (imag(Y));
    %YSignreal = Yreal;
    %YSignimg = Yimg;
    YSign = YSignreal + 1i * YSignimg;
    YSign = Y;
   %xdecoded = zeros(size(x));
   %xdecoded = zeros(nt, size(x,2)*2);
   %symbol = zeros(Nt,2) ;
    
    if strcmp(ModulationType,'BPSK')
        M= 2;
        l = log2(M);
        stored_letters = 0:M-1;
        stored_complex = exp(1i * 2 * pi * stored_letters / M);
            stored_symbols = permn(stored_complex, Nt).';
        
    elseif strcmp(ModulationType,'QPSK')
        M= 4;
        l = log2(M);
        stored_letters = 0:M-1;
        stored_complex = exp(1i * 2 * pi * stored_letters / M);% .* exp(1i *pi/4);
        stored_symbols = permn(stored_complex, Nt).';
    elseif strcmp(ModulationType,'8-PSK')
        M= 8;
        l = log2(M);
        stored_letters = 0:M-1;
        stored_complex = exp(1i * 2 * pi * stored_letters / M);% .* exp(1i *pi/4);
        stored_symbols = permn(stored_complex, Nt).';
        stored_symbols_STF = zeros(Nt, size(stored_symbols,2)*2);
    elseif strcmp(ModulationType,'16-PSK')
        M= 16;
        l = log2(M);
        stored_letters = 0:M-1;
        stored_complex = exp(1i * 2 * pi * stored_letters / M);% .* exp(1i *pi/4);
        stored_symbols = permn(stored_complex, Nt).';
        stored_symbols;
    elseif strcmp(ModulationType,'16-QAM')
        M= 16;
        l = log2(M);
        stored_letters = 0:M-1;
        B = de2bi(stored_letters);
        b = reshape(B', [],1);
        stored_complex = modulate(b, ModulationType );
        stored_symbols = permn(stored_complex, Nt).';

    end
    snr =  SNR * 1;
    %snr =  SNR * sqrt(1 + (1^2));
    snr = sqrt(1)/(sqrt(2)*snr);

    
    fminMod = 0;
    fmin_mod = zeros(size(stored_complex,2),N);
    for u =1:U
        
        rcvd_stat1 = zeros(1,N);
        rcvd_stat2 = zeros(1,N);

        h1 = squeeze(H_est(1,:,u));
        h2 = squeeze(H_est(2,:,u));
       

        rx = YSign(:,u).';
        rcvd_pow = 0;
                
            
            rcvd_stat1(1:2:end) = rcvd_stat1(1:2:end)  + (snr*(conj(h1(1:2:end)).*rx( 1:2:end)...
            +h2(2:2:end).*conj(rx( 2:2:end )))); %for s1
            
        
            rcvd_stat1(2:2:end) = rcvd_stat1(2:2:end)  + (snr*(conj(h2(1:2:end)).*rx( 1:2:end)...
            -h1(2:2:end).*conj(rx(2:2:end)))); %for s2
        
            rcvd_pow = rcvd_pow + (abs(h1).^2+abs(h2).^2);

            
            for ppp = 1:size(stored_complex,2)
                symbol = stored_complex(ppp);
                fmin_mod(ppp,:) = fmin_mod(ppp,:) + (abs(rcvd_stat1 - symbol).^2);
            end
       

    end
    
    %{
    for n = 1:N
        [~,argc] = min((fmin_mod(:,n)));

        shat(n,:) = stored_symbols(:,argc);
    end
    %}
    [~,argc] = min(fmin_mod);
    shat =  stored_complex(:,argc);
        
    decodedSym = shat;

end