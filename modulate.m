
function modulatedSymbols = modulate(x, ModulationType )
%% This function performs modulation.
%% BPSK, 8-PSK 16-PSK modulation types are allowed.
%% 16-QAM and 64-QAM are also allowed.

    
    if nargin <= 1 
        ModulationType = 'BPSK'; %BPSK, 8-PSK, 16-PSK;
    end
    
    if strcmp(ModulationType,'BPSK')
        M= 2;
        modulatedSymbols = MyPSK(x, M).';
    elseif strcmp(ModulationType,'QPSK')
         M= 4;
        modulatedSymbols = MyPSK(x, M).';
    elseif strcmp(ModulationType,'8-PSK')
         M= 8;
        modulatedSymbols = MyPSK(x, M).'; 
    elseif strcmp(ModulationType,'16-PSK')
         M= 16;
        modulatedSymbols = MyPSK(x, M).';
    elseif strcmp(ModulationType,'16-QAM')
        symorder = 'gray';
        M=16;
        modulatedSymbols = qammod(x,M,symorder,'InputType','bit','UnitAveragePower', true,'PlotConstellation',false).';
    elseif strcmp(ModulationType,'64-QAM') 
        symorder = 'gray';
        M=64;
        modulatedSymbols = qammod(x,M,symorder,'InputType','bit','UnitAveragePower', true,'PlotConstellation',true).';
    elseif strcmp(ModulationType,'256-QAM') 
        symorder = 'gray';
        M=256;
        modulatedSymbols = qammod(x,M,symorder,'InputType','bit','UnitAveragePower', true,'PlotConstellation',true).';
    elseif strcmp(ModulationType,'1024-QAM') 
        symorder = 'gray';
        modulatedSymbols = qammod(x,M,symorder,'InputType','bit','UnitAveragePower', true,'PlotConstellation',true).';
    else
        error('Incompatible modulation tyoe');
    end
end