
% differenrial encoder for space time code
% in: this is a cell contains the input STC codewords
% out: Differeniral codes
function [precodedSymbol,NAntennaSym]=spacetimeEncoder(in,Nt,NSub, phasePilots)

    %out(1) = 1;
    in= in;
    N = size(in,2);
    %Nt  = size(in,1);
    NAntennaSym = N/Nt; 
    Ncp =64;
    out = zeros(Nt,1);
    m = 1;
    
    
    N_Spacing = 0;
    N_Carriers = size(in,2) +length(phasePilots);
    N_Pilots = length(phasePilots);
    N_Data = N_Carriers - N_Pilots - 2*N_Spacing;
    N_cp = 32;
    PilotSpacing = ceil((N_Carriers -(2*N_Spacing)) / (N_Pilots));


    pilot_ind = 1+N_Spacing:PilotSpacing:N_Carriers - N_Spacing -1;
    hold_ind = 1:1:N_Carriers;
    Data_Ind  = setdiff(hold_ind, pilot_ind);
    Data_Ind = Data_Ind(1+N_Spacing:length(Data_Ind) -N_Spacing );
    
   
    
    
    precodedSymbol = zeros(Nt, N_Carriers);
    %precodedSymbol(1,1:2:end) = in(1:2:end);
    %precodedSymbol(1,2:2:end) = -conj(in(2:2:end));
    %precodedSymbol(2,1:2:end) = in(2:2:end);
    %precodedSymbol(2,2:2:end) = conj(in(1:2:end));
    
    precodedSymbol(1,Data_Ind) = in(1,:);
    precodedSymbol(2,Data_Ind) = in(2,:);

    precodedSymbol(1,pilot_ind(1:2:end)) = phasePilots(1:2:end);
    precodedSymbol(2,pilot_ind(2:2:end)) = phasePilots(1:2:end);
    precodedSymbol = precodedSymbol./sqrt(2);
    

    %precodedSymbol =  ifft(precodedSymbol,[],2);
    %powerxEncoded = sqrt(var(precodedSymbol,0,2)) ;
    %precodedSymbol = bsxfun(@rdivide,precodedSymbol,powerxEncoded(:)) ;
    %precodedSymbol = precodedSymbol./sqrt(2);


 %precodedSymbol =[precodedSymbol(:, size(precodedSymbol,2) - Ncp+1 :size(precodedSymbol,2) )  precodedSymbol];
 precodedSymbol;
 
 

 
end


