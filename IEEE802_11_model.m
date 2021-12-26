function PDP=IEEE802_11_model(sigma_t,Ts)
    % IEEE 802.11 channel model PDP generator
    % Input:
    % sigma_t : RMS delay spread
    % Ts : Sampling time
    % Output:
    % PDP : Power delay profile
    lmax = ceil(10*sigma_t/Ts); % Eq.(2.6)
    sigma02=(1-exp(-Ts/sigma_t))/(1-exp(-(lmax+1)*Ts/sigma_t)); % Eq.(2.9)
    l=0:lmax; PDP = sigma02*exp(-l*Ts/sigma_t); % Eq.(2.8)
    PDP;
end