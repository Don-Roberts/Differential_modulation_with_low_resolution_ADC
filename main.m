
%% Detects phase information 
%% Equation 40 - ML1BitDetector
%% Equation 53 - ML1BitDetectorRedOFDM


clc
clear
close all

%%
rate =0.5;

Nt = 2;

% CQI = 0.5 set at Nr = 1
Nr = 1;
%CQI = 0.5; 

% CQI = 0.1 set at Nr = 2
%Nr = 2;
CQI = 0;


%Nr = 1;
%CQI = 0.5; 

U = 48; 


%Ts = 1/(20 * 10^6);
%t_rms = 1/(1.5*10^7);

scale = 10^-9;
Ts = 50 *scale;
Ncp = 32;
NcpNoDiff =32;
channelEstType = 1;


fd =5;
%fd =50;
%fd =500;
%fd =5000;
%fd =50000;

pilotFrac = 3;
trmsValue = 50; % 50
%trmsValue = 100; % 50
%trmsValue = 150; % 50

%t_rms = 15 * scale;
%t_rms = 75 * scale;
%t_rms = 100 * scale;
t_rms = trmsValue * scale;

%t_rms = 200 * scale;
%t_rms = 400 * scale;

num_ch = 100; % Number of channels
NSub=64; % FFT size

channeslSteps = 9;
decodingBlocks = 10;
Nb = Nt*100*channeslSteps;               % number of bits per run
Iter = 20;
SNRdB = linspace(-5,10,4);%linspace(0,40,21);
%SNRdB = [15];
SNRList = 10.^(SNRdB/10);
%SNRList = [1];

ModulationType = 'QPSK';
M = 4; 
bitsPerSymbol = log2(M); 
errMode = 'SER';
Nb = Nt  *NSub; 
Nb1 = Nb -2;

ModulationTypeNoDiff = 'QPSK';
MNoDiff = 4; 
bitsPerSymbolNoDiff = log2(MNoDiff); 


NbNoDiff_initial_pilots = 64;
NbNoDiff_initial = 32;

%% pilotFrac

if (pilotFrac == 1)
    NbNoDiff_initial_pilots = 16;
    NbNoDiff_initial = 56;
elseif (pilotFrac == 2)
    NbNoDiff_initial_pilots = 32;
    NbNoDiff_initial = 224; 
elseif (pilotFrac == 3)
    NbNoDiff_initial_pilots = 64;
    NbNoDiff_initial = 192; 
elseif (pilotFrac == 4)
    NbNoDiff_initial_pilots = 128;
    NbNoDiff_initial = 128; 
elseif (pilotFrac == 5)
    %NbNoDiff_initial_pilots = 96;
    %NbNoDiff_initial = 16; 
end

    

x = randi([0 M-1],Nb1,1);
xNoDiff = randi([0 MNoDiff-1],NbNoDiff_initial*2,1);
phaseX = randi([0 M-1],Nb,1);
%x =[0 0 0 1 1 0 1 1];

B = de2bi(x); 
BNoDiff = de2bi(xNoDiff); 
phaseB =  de2bi(phaseX);

b = reshape(B', [],1);
bNoDiff = reshape(BNoDiff', [],1);
bphase = reshape(phaseB', [],1);

S  = modulate(b, ModulationType ).';
SNoDiff  = modulate(bNoDiff, ModulationTypeNoDiff ).';
Sphase = modulate(bphase, ModulationType ).';

%S = 2.*b-1;
ofdmSymbol = 5;
ofdmPilot = 2;
encodedSym  =S;
differentialEncod = diff_encoder(encodedSym,Nt,NSub);
differentialEncod = repmat(differentialEncod,1,ofdmSymbol);
superdifferentialEncod = zeros(size(differentialEncod));
%S_ = reshape(encodedSym , Nt,[]);
%overflow = overflow/Nt;
%Nb_ = Nt*(Nb+overflow);
S_ = SNoDiff.';
symbolsAlamouti = zeros(Nt,NbNoDiff_initial*2);
symbolsAlamouti(1,1:2:end) = S_(1:2:end);
symbolsAlamouti(1,2:2:end) = -conj(S_(2:2:end));
symbolsAlamouti(2,1:2:end) = S_(2:2:end);
symbolsAlamouti(2,2:2:end) = conj(S_(1:2:end));


Sphase = ones(1,NbNoDiff_initial_pilots);
[S_differentialEncod, pilotSymbols]=spacetimeEncoder(symbolsAlamouti,Nt,length(S_),Sphase);
S_differentialEncod = repmat(S_differentialEncod,1,ofdmSymbol);
pilotSymbols = repmat(pilotSymbols,1,ofdmPilot);
    Sphase = Sphase;

%S_differentialEncod = [pilotSymbols S_differentialEncod ];
Nb_ = Nb;
Nb_ = size(differentialEncod,2);
Nb_NoDiff = size(S_differentialEncod,2);

t = 0:Ts: Ts* (Nb_ -1); 
t_NoDiff = 0:Ts: Ts* (Nb_NoDiff -1); 

NsubPaths = 40;
delay =zeros(Ncp ,1);
delay(1) = 1;
lenDela = 0:19;


NbNoADC  = length(S_) +  length(Sphase);
 NtNoADC = 1;
%%
H = zeros(Nt,Nb_, U)/sqrt(2) + 1i/sqrt(2)*zeros(Nt,Nb_, U);
Y = zeros(Nb_, U)/sqrt(2) + 1i/sqrt(2)*zeros(Nb_, U);

YNoADC = zeros(Nb_NoDiff, U)/sqrt(2) + 1i/sqrt(2)*zeros(Nb_NoDiff, U);
H_ = zeros(Nt,Nb_NoDiff, U)/sqrt(2) + 1i/sqrt(2)*zeros(Nt,Nb_NoDiff, U);
H_est = zeros(NtNoADC,Nb_NoDiff, U)/sqrt(2) + 1i/sqrt(2)*zeros(NtNoADC,Nb_NoDiff, U);
H_estReal = zeros(NtNoADC,Nb_NoDiff, U)/sqrt(2) + 1i/sqrt(2)*zeros(NtNoADC,Nb_NoDiff, U);
H_estImag = zeros(NtNoADC,Nb_NoDiff, U)/sqrt(2) + 1i/sqrt(2)*zeros(NtNoADC,Nb_NoDiff, U);


H_estfftsymbol = zeros(2,NbNoADC,U);

%H_est = zeros(Nt,(ofdmPilot + ofdmSymbol)*NbNoADC,U);
%H_estReal = zeros(Nt,(ofdmPilot + ofdmSymbol)*NbNoADC,U);
%H_estImag = zeros(Nt,(ofdmPilot + ofdmSymbol)*NbNoADC,U);

              
%S = (modulate(B, ModulationType )).';
%1/sqrt(2)*randn(Nr_nr2,Nt_nr2) + 1i/sqrt(2)*randn(Nr_nr2,Nt_nr2);
                                                                                                                                               
formatSpec = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d';
%formatSpecRx = "_k_%d_n_%d_tx_%d_rx_%d_codeStyle_%s_txMode_%s_codeRate_%s_csitnone_decoderMethod_%sconstellationDataRx.mat";

str = sprintf(formatSpec,M,Nt,Nr,U,fd,trmsValue);

formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
%formatSpecRx = "_k_%d_n_%d_tx_%d_rx_%d_codeStyle_%s_txMode_%s_codeRate_%s_csitnone_decoderMethod_%sconstellationDataRx.mat";
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);

for k=1:length(SNRdB)
    k
    if (k > 3)
        Iter = 50;
    end

    biterr = 0;
    sererr = 0;
    
    biterr1bit = 0;
    sererr1bit = 0;
    
    biterr1bitRed = 0;
    sererr1bitRed = 0;
    
    biterr1bitRed = 0;
    sererr1bitRed = 0;
    
    biterr1BitNoDiff = 0;
    sererr1BitNoDiff = 0;
    for i=1:Iter

        if (mod( i,10) == 0)
            i
        end
        for Uind = 1:U 
            noise = sqrt(1/(2*   SNRList(k))) * (randn(Nb_, 1) +...
            1i*randn(Nb_, 1));
        
            noiseNoADC =  sqrt(1/(2*SNRList(k))) * (randn(Nb_NoDiff, 1) +...
            1i*randn(Nb_NoDiff, 1));
        
            PDP_1=IEEE802_11_model(t_rms,Ts);
            PDP_2=IEEE802_11_model(t_rms,Ts);

            delay =zeros(length(PDP_1) ,1);
            Llen = 0:length(PDP_1)-1;

            delay(1) = 1;

            delay(2:end) = sort(randi([2,Ncp-10], length(PDP_1)-1,1));
            
            for nt = 1:Nt
            
                for kkk=1:length(PDP_1)
                    h1(:,kkk) = Ray_model(num_ch).*sqrt(PDP_1(kkk));
                    avg_pow_h_1(kkk)= mean(h1(:,kkk).*conj(h1(:,kkk)));

                    h2(:,kkk) = Ray_model(num_ch).*sqrt(PDP_2(kkk));
                    avg_pow_h_2(kkk)= mean(h2(:,kkk).*conj(h2(:,kkk)));
                end

                %hUser_ = randn(Nt,(Nb)/decodingBlocks)/sqrt(2) + 1i/sqrt(2)*randn(Nt,(Nb)/decodingBlocks);
                hUser1 = fft(h1(1,:),Nb_);
                hUser2 = fft(h2(1,:),Nb_NoDiff);

                H(nt,:, Uind) =  hUser1;%./ sqrt(1*var(hUser1)) ;;
                H_(nt,:, Uind) =  hUser2;%./ sqrt(1*var(hUser2)) ;;
            end 
            
            %TxChan = zeros(size(differentialEncod));
            %TxChanNoDiff = zeros(size(S_differentialEncod));

            %TxChan(:,:) = (PDP_1(1)) * squeeze(channHdim(:,:));
            %YHold = sum(TxChan .* differentialEncod,1);
 
            %TxChanNoDiff(:,:) = (PDP_1(1)) *  squeeze(channHdimNoDiff(:,:));
            %YHold_ =  sum(TxChanNoDiff .* S_differentialEncod,1);
            


            if Nt ~= 1
                YHold = sum(squeeze(H(:,:,Uind)) .*differentialEncod,1).';
                %YHold = sum( differentialEncod,1).';
                %YHold_ = sum(squeeze(H_(:,:,Uind)) .*S_,1).';
            else
                YHold = (squeeze(H(:,:,Uind)) .*differentialEncod).';
                %YHold_ = (squeeze(H_(:,:,Uind)) .*S_).';
            end
            
            if Nt ~= 1
                %YHold_ = sum(squeeze(H_(:,:,Uind)) .*S_,1).';
                YHold_ = sum(squeeze(H_(:,:,Uind)) .*S_differentialEncod,1).';
                %YHold_ = sum(S_differentialEncod,1).';

            else
                %YHold = (squeeze(H(:,:,Uind)) .*differentialEncod).';
                YHold_ = (squeeze(H_(:,:,Uind)) .*S_differentialEncod).';
            end
            

        
            Y(:, Uind) = YHold+noise ;
            YNoADC(:, Uind) = YHold_+ noiseNoADC; 

 
            
            snrEst = sqrt(var(differentialEncod) / var(noise));

            snrEstNoADC = sqrt(var(S_) / var(noiseNoADC));

            snrEstNoADC;
            
                 
                YSignreal = sign(real(Y(:, Uind)));
                YSignimg = sign(imag(Y(:, Uind)));
                YSign(:,Uind) = YSignreal + 1i * YSignimg;
                 %YSign = Y;
                 
                YSignrealNoADC = sign(real(YNoADC(:, Uind)));
                YSignimgNoADC = sign(imag(YNoADC(:, Uind)));
                YSignNoADC(:,Uind) = YSignrealNoADC + 1i * YSignimgNoADC;
                %YSignNoADC = YNoADC;
                
                 YSignShape(:,:, Uind) = reshape(YSign(:,Uind),[],ofdmSymbol).';
                 YShape(:,:, Uind) = reshape(Y(:,Uind),[],ofdmSymbol).';
                 YShapeNoADC(:,:, Uind) = reshape(YSignNoADC(:,Uind),[],ofdmSymbol).';

                 YShape;
                 %YSignShape_cp(:,:, Uind) = YSignShape(:,Ncp+1:size(YSignShape,2),Uind);
                 %YShape_cp(:,:, Uind) = YShape(:,Ncp+1:size(YShape,2),Uind);
                 %YShapeNoADC_cp(:,:, Uind) = YShapeNoADC(:,NcpNoDiff+1:size(YShape,2),Uind);

                 %YSignfft(:,:, Uind) = fft(YSignShape_cp(:,:, Uind) ,[],2);
                 %Yfft(:,:, Uind) = fft(YShape_cp(:,:, Uind),[],2);
                 %YfftNoADC(:,:, Uind) = fft(YShapeNoADC_cp(:,:, Uind),[],2);
               

                 YSignfft(:,:, Uind) = YSignShape(:,:, Uind);
                 Yfft(:,:, Uind) = YShape(:,:, Uind);
                 YfftNoADC(:,:, Uind) = YShapeNoADC(:,:, Uind);
                 
                 %pilotfftNoADC(:,:, Uind) = YfftNoADC(1:ofdmPilot,:,Uind);
                 %YsymbolfftNoADC(:,:, Uind) =  YfftNoADC(ofdmPilot+1:ofdmPilot+ofdmSymbol,:,Uind);
                 
                 %lenDela = 1:length(PDP_1) -4;
             
             
                 N_Spacing = 0;
                N_Pilots = length(Sphase(1:NbNoDiff_initial_pilots));
                N_Data = NbNoADC - N_Pilots - 2*N_Spacing;
                PilotSpacing = (NbNoADC -(2*N_Spacing)) / (N_Pilots);
                PilotSpacing = ((NbNoADC -(2*N_Spacing)) / (N_Pilots));


                pilot_ind = 1+N_Spacing:PilotSpacing:NbNoADC - N_Spacing -1;
                hold_ind = 1:1:NbNoADC;
                Data_Ind  = setdiff(hold_ind, pilot_ind);
                Data_Ind = Data_Ind(1+N_Spacing:length(Data_Ind) -N_Spacing );

             H_estPilotsTemp = zeros(Nt,(ofdmSymbol),length(Data_Ind));

    
            for nt =1:2
 
            
             
             Dnl = (1/sqrt(NbNoADC) ) * exp(-1i * 2 * pi  * [0:NbNoADC-1]'*lenDela/ NbNoADC);
             intex_nt = pilot_ind(nt:2:end);
             for nofdm = 1:(ofdmSymbol)
                 
                 pilotfftNoADC = YfftNoADC(nofdm,intex_nt,Uind);
            %Dnl = (1/sqrt(Nb) ) * exp(-1i * 2 * pi  * [0:Nb-1]'*Llen/ Nb);
            %Q = exp(-1i * 2 * pi * intex'*Llen/ Nb) / sqrt(Nb);
        
            %hLS = inv(Q' * Q) * Q' * pilotfftNoADC(nofdm,nt:2*U:NbNoADC,Uind).';
            %HLs = Dnl * hLS;
            %HLs;
                Q = exp(-1i * 2 * pi * intex_nt'*lenDela/ NbNoADC) / sqrt(NbNoADC);

                hLS = inv(Q' * Q) * Q' * pinv(diag(Sphase(nt:2:end)))*pilotfftNoADC.';
                HLs = Dnl * hLS;
                %intex = [intex ((nt:2:NbNoADC)+((NbNoADC *(nofdm - 1))))];
                H_estPilotsTemp(nt,nofdm,:) = HLs(Data_Ind);
             end
             
         pilotIndY =1:ofdmSymbol;

        dataIndXPlane = repmat(Data_Ind, ofdmSymbol ,1);
        dataIndYPlane = repmat(1:ofdmSymbol,length(Data_Ind),1)';

        pilotIndXPlane = repmat(Data_Ind, ofdmSymbol,1);
        pilotIndYPlane = repmat(pilotIndY, length(Data_Ind),1)';

        hData(nt,:,:,Uind) = squeeze(H_estPilotsTemp(nt,:,:));
        %hData(nt,:,:,Uind) = interp2(pilotIndXPlane,pilotIndYPlane,squeeze(H_estPilotsTemp(nt,:,:)),dataIndXPlane,dataIndYPlane,'spline');;
        hData;
        end
    Y;
       %% Plots
         
        end
             %H_est = H_estReal + (1i*H_estImag) ;
            %H_estFinal = H_est(:,ofdmPilot*NbNoADC+1:(ofdmPilot+ofdmSymbol)*NbNoADC,:);
     for nofdm = 1:ofdmSymbol
     
         Yfftsymbol = squeeze(Yfft(nofdm,:,:));
         YSignfftsymbol = squeeze(YSignfft(nofdm,:,:));
         YNoADCfftsymbol = squeeze(YfftNoADC(nofdm,Data_Ind,:));
         %H_estfftsymbol = squeeze(H_estfft(nofdm,:,:,:));
         %H_estfftsymbol = H_est(:,((nofdm-1)*NbNoADC)+1: (nofdm)*NbNoADC,:);
         %H_estfft;
         %H_estfftsymbol = H_estFinal(:,((nofdm-1)*NbNoADC)+1: (nofdm)*NbNoADC,:);
         H_estfftsymbol = squeeze(hData(:, nofdm,:,:)) ;
         H_estfftsymbol;
                
     
     %shat1bit = ML1BitDetectorOFDM(YSignfftsymbol ,Nb, Nt,ModulationType,snrEst);
     %bhat1bit = demodulate(shat1bit, ModulationType );
     %biterr1bit = biterr1bit+sum(bhat1bit ~= b);
     %biterr1bit = 0;
     
     %shat1bit = bi2de(reshape(bhat1bit, bitsPerSymbol ,[]).');
     %sererr1bit = sererr1bit+sum(shat1bit ~= x);
     %sererr1bit = 0;
     
      if (k <4)
     snrEst = sqrt(var(noiseNoADC)/ (1*var(S_differentialEncod(1,:))) );
     current = snrEst;
    else
        snrEst = current;
     snrEst = sqrt(var(noiseNoADC)/ (1*var(S_differentialEncod(1,:))) );

      end
    
      % One bit differential detector
      shat1bit = ML1BitDetectorOFDM(YSignfftsymbol ,Nb, Nt,ModulationType,snrEst);
     bhat1bit = demodulate(shat1bit, ModulationType );
     biterr1bit = biterr1bit+sum(bhat1bit ~= b);
     
     shat1bit = bi2de(reshape(bhat1bit, bitsPerSymbol ,[]).');
     sererr1bit = sererr1bit+sum(shat1bit ~= x);
     sererr1bit;
     
     %One bit differential detector with reduced complexity
      shat1bitRed = ML1BitDetectorRedOFDM(YSignfftsymbol ,Nb, Nt,ModulationType,snrEst);
     bhat1bitRed = demodulate(shat1bitRed, ModulationType );
     biterr1bitRed = biterr1bitRed+sum(bhat1bitRed ~= b);
     
     shat1bitRed = bi2de(reshape(bhat1bitRed, bitsPerSymbol ,[]).');
     sererr1bitRed = sererr1bitRed+sum(shat1bitRed ~= x);
     sererr1bitRed;
     
     
     % Coherent Alamouti Scheme With One Bit
     shat1BitNoDiff = ML1BitAlamouti(YNoADCfftsymbol, NbNoADC, NtNoADC,H_estfftsymbol,ModulationTypeNoDiff,snrEst);
     bhat1BitNoDiff = demodulate(shat1BitNoDiff, ModulationTypeNoDiff );
     biterr1BitNoDiff = biterr1BitNoDiff+sum(bhat1BitNoDiff ~= bNoDiff);
     
     shat1BitNoDiff = bi2de(reshape(bhat1BitNoDiff, bitsPerSymbol ,[]).');
     sererr1BitNoDiff = sererr1BitNoDiff+sum(shat1BitNoDiff ~= xNoDiff);
     sererr1BitNoDiff;
          
     
     shat = MLDetector(Yfftsymbol ,Nb, Nt,ModulationType,snrEst);
     bhat = demodulate(shat, ModulationType );
     biterr = biterr+sum(bhat ~= b);
     
     shat = bi2de(reshape(bhat, bitsPerSymbol ,[]).');
     sererr = sererr+sum(shat ~= x);
     sererr;
     
     %{
     shat = MLDetector(Yfftsymbol ,Nb, Nt,ModulationType,snrEst);
     bhat = demodulate(shat, ModulationType );
     biterr = biterr+sum(bhat ~= b);
     
     shat = bi2de(reshape(bhat, bitsPerSymbol ,[]).');
     sererr = sererr+sum(shat ~= x);
     sererr;
     %UnAttshat1bit = UnAttML1BitDetector(Y ,Nb, Nt,ModulationType,snrEst);
     %UnAttbhat1bit = demodulate(UnAttshat1bit, ModulationType );
     %UnAttbiterr1bit = UnAttbiterr1bit+sum(UnAttbhat1bit ~= b);
     
     %UnAttshat1bit = bi2de(reshape(UnAttbhat1bit, bitsPerSymbol ,[]).');
     %UnAttsererr1bit = UnAttsererr1bit+sum(UnAttshat1bit ~= x);
     %UnAttsererr1bit;
     
     
     shatRed = MLDetectorRed(Yfftsymbol ,Nb, Nt,ModulationType,snrEst);
     bhatRed = demodulate(shatRed, ModulationType );
     biterrRed = biterrRed+sum(bhatRed ~= b);
     
     shatRed = bi2de(reshape(bhatRed, bitsPerSymbol ,[]).');
     sererrRed = sererrRed+sum(shatRed ~= x);
     sererrRed;
     %}

     %shat1BitNoDiff = ML1Bit(Y ,Nb, Nt,ModulationType,snrEst);
     %H_est = H_ + sqrt(CQI) *( randn(size(H_)) + (1i*randn(size(H_))));
     
     %H_est = H_ + sqrt(CQI) *( randn(size(H_)) + (1i*randn(size(H_))));

     %{
     shat1BitNoDiff = ML1Bit(YNoADC, NbNoADC, NtNoADC,H_est,ModulationType,snrEst);
     bhat1BitNoDiff = demodulate(shat1BitNoDiff, ModulationType );
     biterr1BitNoDiff = biterr1BitNoDiff+sum(bhat1BitNoDiff ~= b);
     
     shat1BitNoDiff = bi2de(reshape(bhat1BitNoDiff, bitsPerSymbol ,[]).');
     sererr1BitNoDiff = sererr1BitNoDiff+sum(shat1BitNoDiff ~= x);
     sererr1BitNoDiff;
     %}
     end
     
   

    end

    
    ber1bit(k) = biterr1bit/length(b)/Iter/ofdmSymbol;
    ber1bitRed(k) = biterr1bitRed/length(b)/Iter/ofdmSymbol
    ber1Diff(k) = biterr/length(b)/Iter/ofdmSymbol
    ber1NoADC(k) = biterr1BitNoDiff/length(bNoDiff)/Iter/ofdmSymbol

    ser1bit(k) = sererr1bit/length(x)/Iter/ofdmSymbol;
    ser1bitRed(k) = sererr1bitRed/length(x)/Iter/ofdmSymbol
    ser1Diff(k) = sererr/length(x)/Iter/ofdmSymbol
    ser1NoADC(k) = sererr1BitNoDiff/length(xNoDiff)/Iter/ofdmSymbol
    
    %perfectber1bit(k) = perfectbiterr1bit/length(b)/Iter
    %newber1bitZF(k) = newZFbiterr1bit/length(b)/Iter
    %ber1bitDiff(k) = biterr1bitDiff/length(b)/Iter   

    %ser1bit(k) = sererr1bit/Iter/length(x);
    %perfectser1bit(k) = perfectsererr1bit/Iter/length(x);
    %newser1bitZF(k) = newZFsiterr1bit/Iter/length(x);

    %ser1bitDiff(k) = sererr1bitDiff/length(x)/Iter  ;  
    
end


%%

figure;
semilogy(SNRdB, ber1Diff,'r-p')
semilogy(SNRdB, ber1bitRed,'g-p')

hold on
grid on; 
semilogy(SNRdB, ber1NoADC,'b-p')
semilogy(SNRdB, ber1bit,'k-p')
%semilogy(SNRdB, UnAttber1bit,'r-p')
%semilogy(SNRdB, ber1bitRed,'g-p')
%semilogy(SNRdB, ber1NoADC,'m-p')

%semilogy(SNRdB, ber1bitDiff,'r-p')

hold off;

xlabel('EbNo (dB)' , 'FontSize', 20, 'Interpreter', 'latex')
ylabel('BER', 'FontSize', 20, 'Interpreter', 'latex')
%legend({'OSTBC : BPSK'; 'NN Based Encoder: QPSK'}, 'Location', 'northeast', 'FontSize', 15, 'Interpreter', 'latex')
%title('Rate: 1 bit/s/Hz  ', 'FontSize', 12, 'Interpreter', 'latex')

%%

%save(strcat('DataOneTime_109_25/ADC/' ,str, '.mat'),'ber1bit')
%save(strcat('Data/ADC_Perfect/' ,str, '.mat'),'perfectber1bit')
%save(strcat('Data/ADC_ZF/' ,str, '.mat'),'newber1bitZF')

%save(strcat('DataOneTime_109_25/FullResolution/' ,str, '.mat'),'ber1')
%save(strcat('DataOneTime_109_25/UnAttADC/' ,str, '.mat'),'UnAttber1bit')
save(strcat('Data/ADC/' ,str, '.mat'),'ber1bit')
save(strcat('Data/ADCRed/' ,str, '.mat'),'ber1bitRed')
save(strcat('Data/FullResolutionDiffSer/' ,str, '.mat'),'ber1Diff')
save(strcat('Data/NoADC/' ,strNoADC, '.mat'),'ber1NoADC')

save(strcat('Data/ADCSer/' ,str, '.mat'),'ser1bit')
save(strcat('Data/ADCRedSer/' ,str, '.mat'),'ser1bitRed')
save(strcat('Data/FullResolutionDiffSer/' ,str, '.mat'),'ser1Diff')
save(strcat('Data/NoADCSer/' ,strNoADC, '.mat'),'ser1NoADC')
