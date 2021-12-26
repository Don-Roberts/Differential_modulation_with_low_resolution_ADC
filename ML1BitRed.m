
function decodedSym=ML1Bit(Y, Nb, Nt,H_,ModulationType,SNR)
    %-1 = 1;
    % 1 = 0
    T = 10;
    N = size(Y,1);
    U = size(Y,2);
    H = eye(2);
    
    %H_ = sign(H_);
    channelSteps = 10;
    
    shat = zeros(Nb,Nt);
    m= 1;

    Yreal = (real(Y));
    Yimg = (imag(Y));
    
    YSignreal = sign(real(Y));
    YSignimg = sign(imag(Y));
    
    
    %ySign = [YSignreal ; YSignimg];
    %YSignreal = Yreal;
    %YSignimg = Yimg;
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
    elseif strcmp(ModulationType,'16-QAM')
        M= 16;
        l = log2(M);
        stored_letters = 0:M-1;
        B = de2bi(stored_letters);
        b = reshape(B', [],1);
        stored_complex = modulate(b, ModulationType );
        stored_symbols = permn(stored_complex, Nt).';
    end
    
    snr = SNR;% sqrt((l*SNR));
    %snr=sqrt(2);
    snr = (1*SNR)/5;% sqrt((l*SNR));
    snr = sqrt(2)*(1*SNR);
    snr = 1;
    for n = 1:N
        fminMod = 0;
        for u =1:U
            fminReal_ = (real(stored_symbols));% .* YSignreal(n-1,u);
            fminImg_ = (imag(stored_symbols));% .* YSignimg(n-1,u);
            %fminImg_ = [0 0];
            %fmin_ = fminReal_ + 1i*fminImg_;
            fmin_ = [fminReal_; fminImg_];
            Hinst = [real(H_(:,n,u))  imag(H_(:,n,u)); -imag(H_(:,n,u)) real(H_(:,n,u))].';
            Hinst_(1,:) = YSignreal(n,u) * Hinst(1,:);
            Hinst_(2,:) = YSignimg(n,u) * Hinst(2,:);

            %newY1 = Hinst_(1,:) * fmin_.';
            %newY2 = Hinst_(2,:)* fmin_.';
            newY = (Hinst_) * fmin_;
            newY;
            %Hinst = [real(H_(:,n,u)).'  -imag(H_(:,n,u)).'; imag(H_(:,n,u)).' real(H_(:,n,u)).'];

            %realY = fminReal_ * YSignreal(n,u);
            %imgY = fminImg_ * YSignimg(n,u);
            %newY1 =  YSignreal(n,u) *fminReal_ - YSignimg(n,u) *fminImg_;
            %newY2 = YSignimg(n,u) *fminReal_   + YSignreal(n,u) * fminImg_;
                     


            fminMod = fminMod + log(0.5*erfc( (-1/sqrt(2))*snr*newY(1,:) ))+...
                log(0.5*erfc( (-1/sqrt(2))*snr*newY(2,:)));
        end
        [~,argc] = max((fminMod));

        shat(m,:) = stored_symbols(:,argc);
        m = m+1;

    end

    %shat = shat -1;
    %shat = 2.*shat -1;
decodedSym  =  reshape(shat.', [],1);
decodedSym;
end