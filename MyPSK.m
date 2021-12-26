function ModulatedSymbols = MyPSK(x, M)
%%Don-Roberts Emenonye

%Input - vector of bits, M
%Output - Vector of complex numbers
%Use Gray coding and unit energy
 
%x = randi([0,1],720,1);
%M = 16;
l = log2(M);

if mod( length(x), l) ~= 0
    p = ceil(length(x) /l);
    pp =  zeros(1,p * l - length(x));
    x = [x  pp] ;
end


len = length(x);

%b_y = reshape(x,len/l,l);
b_y = reshape(x,l,len/l)';

b_hold(:,1) = b_y(:,1);
for i = 2:l
    b_hold(:,i) = xor(b_y(:,i-1), b_y(:,i));
end
b_y = b_hold;



factors_2 = 2.^([l-1:-1:0]);

dft = zeros(length(x)/l, 1);

for i=1:l
    dft = dft + factors_2(i)*b_y(:,i);
end

%ModulatedSymbols = exp(1i*2*pi*dft/M);% .*  exp(1i*(l-1)*pi)  ;
ModulatedSymbols = exp(1i*2*pi*dft/M) ;%.*  exp(1i*pi/8)  ;

y = ModulatedSymbols;

