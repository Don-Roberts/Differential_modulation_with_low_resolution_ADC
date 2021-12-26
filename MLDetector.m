function [decodedSym,shat_mag]=MLDetector(Y,NSub, Nt,ModulationType, SNR)
    %-1 = 1;
    % 1 = 0
    N = size(Y,1);
    U = size(Y,2);

    NSub = N/2;
    shat = zeros(1,2);
    shat_mag = zeros(NSub,1);
    m= 1;

    
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
        stored_symbols;
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
    snr = sqrt((2*l*SNR));
    snr = sym(sqrt(2)*SNR);% sqrt((l*SNR));
    snr = (1*SNR)/5;% sqrt((l*SNR));
    snr = sqrt(2)*(1*SNR);
    snr =1;
    
    
        for p = 1:NSub-1
            fminMod = 0;
            fminMod_1 = 0;
            if  (p == 100000000)
            p;
            else
                
                argHold = -999999999;
                argcHold = -99;
                fmin_mod = zeros(1, size(stored_symbols,2));
                %fminMod_1 = zeros(3, length(stored_symbols_mag));
                
                %Yinst1 = zeros(2,Nt);
                %Yinst2 = zeros(2,Nt);
                
                %YBlock1 = Y((p-2)*4+1:(p-1)*4,:);
                %YBlock2 = Y((p-1)*4+1:(p)*4,:);
                
                Yinst1 = zeros(2,1);
                Yinst2 = zeros(2,1);
                
                YBlock1 = Y(2*p-1 :2*p,:);
                YBlock2 = Y(2*(p+1)-1 :2*(p+1),:);
                %2*p-1 :2*p
                %2*(p+1)-1 :2*(p+1)
                for u =1:U
                    
                    %Yinst1(:,1) = YBlock1(1:2:4,u);
                    %Yinst1(:,2) = YBlock1(2:2:4,u);

                    %Yinst2(:,1) = YBlock2(1:2:4,u);
                    %Yinst2(:,2) = YBlock2(2:2:4,u);
                    
                    Yinst1 = YBlock1(:,u);
                    Yinst2 = YBlock2(:,u);

                    for ppp = 1:size(stored_symbols,2)
                        symbol = zeros(Nt,2) ;
                        symbol(:,1) = stored_symbols(1:Nt,ppp);
                        symbol(:,2) = conj(stored_symbols(flip(1:2),ppp));
                        symbol(1,2) = -1*symbol(1,2);
                        symbol =   (1/sqrt(2))*symbol;
                        fmin_mod(ppp) = fmin_mod(ppp) + sum(abs(Yinst2 - (symbol.' * Yinst1)).^2);

                    end
                    
                    

                end
            [~,argc] = min((fmin_mod));

            
            %    shat_mag(m) = ;
            %argc=1;
            shat(m,:) = stored_symbols(:,argc);
            m = m+1;
            end
        end

    %shat = shat -1;
    %shat = 2.*shat -1;
decodedSym  =  reshape(shat.', [],1);
%decodedSym = shat;
decodedSym =  reshape(shat, [],1);;
%decodedSym = shat;

end