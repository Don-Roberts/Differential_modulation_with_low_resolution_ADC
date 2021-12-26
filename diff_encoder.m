% differenrial encoder for space time code
% in: this is a cell contains the input STC codewords
% out: Differeniral codes
function [differentialEncod]=diff_encoder(in,Nt,NSub)

    %out(1) = 1;
    in= in.';
    N = size(in,2);
    %Nt  = size(in,1);
    NAntennaSym = N/Nt; 
    Ncp =32;
    out = zeros(Nt,1);
    m = 1;
    
    x = zeros(Nt,NAntennaSym);
    for nt = 1:Nt
        x( nt,:) = in(((nt-1)*NAntennaSym) +1:(nt)*NAntennaSym);
    end
    
   xEncoded = zeros(size(x));
   xEncoded = zeros(nt, size(x,2)*2);

   symbol = zeros(Nt,2) ;
    
    for p = 1:NSub
       symbol = zeros(Nt,2) ;

        if (p == 1)
           symbol = eye(Nt,2) ;
            xEncoded(:,(p-1)*2+1:(p)*2) = (1/1)*symbol;

        else
           symbol(:,1) = x(1:Nt,p-1);
           symbol(:,2) = conj(x(flip(1:2),p-1));
           symbol(1,2) = -1*symbol(1,2);
            xEncoded(:,(p-1)*2+1:(p)*2) = xEncoded(:,(p-2)*2+1:(p-1)*2)* (1/sqrt(2))*symbol;
            %(p-1)*2+1:(p)*2
        end
        xEncoded;
    end
    
    
    xEncoded;
    
 differentialEncod=zeros(size(xEncoded));
 differentialEncod = xEncoded;
% differentialEncod(1,1:NSub) = xEncoded(1,1:2:NSub);
% differentialEncod(2,NSub+1:2*Nsub) = xEncoded(2,2:2:NSub);

%xEncoded =  ifft(xEncoded,[],2);
%powerxEncoded = sqrt(var(xEncoded,0,2)) ;
%xEncoded = bsxfun(@rdivide,xEncoded,powerxEncoded(:)) ;
%xEncoded = xEncoded./sqrt(2);
 %differentialEncod =[xEncoded(:, size(xEncoded,2) - Ncp+1 :size(xEncoded,2) )  xEncoded];
 differentialEncod;
end