function deModulatedSymb = demodulate(symbols, deModulationType )
%% This function performs modulation.
%% BPSK, 8-PSK 16-PSK modulation types are allowed.
%% 16-QAM and 64-QAM are also allowed.

    
    if nargin <= 1 
        deModulationType = 'BPSK'; %BPSK, 8-PSK, 16-PSK;
    end
    
    
    if strcmp(deModulationType,'BPSK')
        M= 2;
        z=MyDetectPSK(symbols,M);
    elseif strcmp(deModulationType,'QPSK')
         M= 4;
        z=MyDetectPSK(symbols,M);
    elseif strcmp(deModulationType,'8-PSK')
         M= 8;
        z=MyDetectPSK(symbols,M); 
        
    elseif strcmp(deModulationType,'16-PSK')
         M= 16;
        z = MyDetectPSK(symbols, M);
        
    elseif strcmp(deModulationType,'16-QAM')
        M= 16;
        z = qamdemod(symbols,M,'gray','UnitAveragePower', true,'OutputType','bit');
        %deModulatedSymbols=MyDetectQAM(symbols,M);         
    elseif strcmp(deModulationType,'64-QAM') 
        M= 64;
        z = qamdemod(symbols,M,'gray','UnitAveragePower', true,'OutputType','bit');
        %deModulatedSymbols=MyDetectQAM(symbols,M);
    elseif strcmp(deModulationType,'256-QAM') 
         M= 256;
         z = qamdemod(symbols,M,'gray','UnitAveragePower', true,'OutputType','bit');
        %deModulatedSymbols=MyDetectQAM(symbols,M);
    elseif strcmp(deModulationType,'1024-QAM') 
         M= 1024;
         z = qamdemod(symbols,M,'gray','UnitAveragePower', true,'OutputType','bit');
        %deModulatedSymbols=MyDetectQAM(symbols,M);
    else
        error('Incompatible de-modulation tyoe');
    end    
    deModulatedSymb =  reshape(z, [],1);
end