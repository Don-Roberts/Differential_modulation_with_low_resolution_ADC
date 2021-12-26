clc
clear
close all

ser_th = 0.05;
bandwidth = 20*10^6;
bandwidth = 1;

number_of_OFDM = 5;
Number_of_subcarriers = 256;

Nt = 2;
Nr = 1;
SNRdB = linspace(0,15,4);%linspace(0,40,21);
SNRdB = linspace(-5,15,5);%linspace(0,40,21);

SNRList = 10.^(SNRdB/10);


ModulationType = '8-PSK';
M = 8; 
bitsPerSymbol = log2(M);  
pilotUsedChannelEst = (256 - 2) * bitsPerSymbol;

fd = 5;
trmsValue = 50;
U = 48;
alowhigh = 2;
formatSpec = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d';
str = sprintf(formatSpec,M,Nt,Nr,U,fd,trmsValue);
ber1bitU128_ = load("Data/ADCRedSer/" +str+ ".mat");
ber1bitRed_48_50 = ber1bitU128_.ser1bitRed;
ber1bitRed_48_50(ber1bitRed_48_50 > ser_th) = 1;
Number_of_subcarriers_used_no_diff = Number_of_subcarriers;
effNoDiff = pilotUsedChannelEst/Number_of_subcarriers;
throughputNoDiff_48_50 = (effNoDiff) *(1-ber1bitRed_48_50) * bandwidth;

U = 64;
alowhigh = 2;
formatSpec = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d';
str = sprintf(formatSpec,M,Nt,Nr,U,fd,trmsValue);
ber1bitU128_ = load("Data/ADCRedSer/" +str+ ".mat");
ber1bitRed_64_50 = ber1bitU128_.ser1bitRed;
ber1bitRed_64_50(ber1bitRed_64_50 > ser_th) = 1;
effNoDiff = pilotUsedChannelEst/Number_of_subcarriers;
throughputNoDiff_64_50 = (effNoDiff) *(1-ber1bitRed_64_50) * bandwidth;

trmsValue = 100;
U = 48;
alowhigh = 2;
formatSpec = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d';
str = sprintf(formatSpec,M,Nt,Nr,U,fd,trmsValue);
ber1bitU128_ = load("Data/ADCRedSer/" +str+ ".mat");
ber1bitRed_48_100 = ber1bitU128_.ser1bitRed;
ber1bitRed_48_100(ber1bitRed_48_100 > ser_th) = 1;
effNoDiff = pilotUsedChannelEst/Number_of_subcarriers;
throughputNoDiff_48_100 = (effNoDiff) *(1-ber1bitRed_48_100) * bandwidth;



U = 64;
alowhigh = 2;
formatSpec = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d';
str = sprintf(formatSpec,M,Nt,Nr,U,fd,trmsValue);
ber1bitU128_ = load("Data/ADCRedSer/" +str+ ".mat");
ber1bitRed_64_100 = ber1bitU128_.ser1bitRed;
ber1bitRed_64_100(ber1bitRed_64_100 > ser_th) = 1;
effNoDiff = pilotUsedChannelEst/Number_of_subcarriers;
throughputNoDiff_64_100 = (effNoDiff) *(1-ber1bitRed_64_100) * bandwidth;


trmsValue = 150;
U = 48;
alowhigh = 2;
formatSpec = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d';
str = sprintf(formatSpec,M,Nt,Nr,U,fd,trmsValue);
ber1bitU128_ = load("Data/ADCRedSer/" +str+ ".mat");
ber1bitRed_48_150 = ber1bitU128_.ser1bitRed;
ber1bitRed_48_150(ber1bitRed_48_150 > ser_th) = 1;
effNoDiff = pilotUsedChannelEst/Number_of_subcarriers;
throughputNoDiff_48_150 = (effNoDiff) *(1-ber1bitRed_48_150) * bandwidth;



U = 64;
alowhigh = 2;
formatSpec = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d';
str = sprintf(formatSpec,M,Nt,Nr,U,fd,trmsValue);
ber1bitU128_ = load("Data/ADCRedSer/" +str+ ".mat");
ber1bitRed_64_150 = ber1bitU128_.ser1bitRed;
ber1bitRed_64_150(ber1bitRed_64_150 > ser_th) = 1;
effNoDiff = pilotUsedChannelEst/Number_of_subcarriers;
throughputNoDiff_64_150 = (effNoDiff) *(1-ber1bitRed_64_150) * bandwidth;


pilotFrac = 4;
pilotUsedChannelEst = 256 - 32;
pilotUsedChannelEst = 256 - 64;
pilotUsedChannelEst = (256 - 128) * bitsPerSymbol;

NtNoADC =1;
fd = 5;
trmsValue = 50;
U = 48;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADCSer/" +strNoADC+ ".mat");
ber1NoADC_48_50 = ber1bitU128_.ser1NoADC;
ber1NoADC_48_50(ber1NoADC_48_50 > ser_th) = 1;
effNoDiff = pilotUsedChannelEst/Number_of_subcarriers;
throughputNoADC_48_50 = (effNoDiff) *(1-ber1NoADC_48_50) * bandwidth;


U = 64;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADCSer/" +strNoADC+ ".mat");
ber1NoADC_64_50 = ber1bitU128_.ser1NoADC;
ber1NoADC_64_50(ber1NoADC_64_50 > ser_th) = 1;
effNoDiff = pilotUsedChannelEst/Number_of_subcarriers;
throughputNoADC_64_50 = (effNoDiff) *(1-ber1NoADC_64_50) * bandwidth;

trmsValue = 100;
U = 48;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADCSer/" +strNoADC+ ".mat");
ber1NoADC_48_100 = ber1bitU128_.ser1NoADC;
ber1NoADC_48_100(ber1NoADC_48_100 > ser_th) = 1;
effNoDiff = pilotUsedChannelEst/Number_of_subcarriers;
throughputNoADC_48_100 = (effNoDiff) *(1-ber1NoADC_48_100) * bandwidth;


U = 64;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADCSer/" +strNoADC+ ".mat");
ber1NoADC_64_100 = ber1bitU128_.ser1NoADC;
ber1NoADC_64_100(ber1NoADC_64_100 > ser_th) = 1;
effNoDiff = pilotUsedChannelEst/Number_of_subcarriers;
throughputNoADC_64_100 = (effNoDiff) *(1-ber1NoADC_64_100) * bandwidth;



trmsValue = 150;
U = 48;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADCSer/" +strNoADC+ ".mat");
ber1NoADC_48_150 = ber1bitU128_.ser1NoADC;
ber1NoADC_48_150(ber1NoADC_48_150 > ser_th) = 1;
effNoDiff = pilotUsedChannelEst/Number_of_subcarriers;
throughputNoADC_48_150 = (effNoDiff) *(1-ber1NoADC_48_150) * bandwidth;



U = 64;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADCSer/" +strNoADC+ ".mat");
ber1NoADC_64_150 = ber1bitU128_.ser1NoADC;
ber1NoADC_64_150(ber1NoADC_64_150 > ser_th) = 1;
effNoDiff = pilotUsedChannelEst/Number_of_subcarriers;
throughputNoADC_64_150 = (effNoDiff) *(1-ber1NoADC_64_150) * bandwidth;


pilotFrac = 2;
pilotUsedChannelEst = (256 - 32) * bitsPerSymbol;
%pilotUsedChannelEst = 256 - 64;
%pilotUsedChannelEst = 256 - 128;

NtNoADC =1;
fd = 5;
trmsValue = 50;
U = 48;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADCSer/" +strNoADC+ ".mat");
ber1NoADC_48_50 = ber1bitU128_.ser1NoADC;
ber1NoADC_48_50(ber1NoADC_48_50 > ser_th) = 1;
effNoDiff = pilotUsedChannelEst/Number_of_subcarriers;
throughputNoADC_48_50_2 = (effNoDiff) *(1-ber1NoADC_48_50) * bandwidth;


U = 64;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADCSer/" +strNoADC+ ".mat");
ber1NoADC_64_50 = ber1bitU128_.ser1NoADC;
ber1NoADC_64_50(ber1NoADC_64_50 > ser_th) = 1;
effNoDiff = pilotUsedChannelEst/Number_of_subcarriers;
throughputNoADC_64_50_2 = (effNoDiff) *(1-ber1NoADC_64_50) * bandwidth;

trmsValue = 100;
U = 48;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADCSer/" +strNoADC+ ".mat");
ber1NoADC_48_100 = ber1bitU128_.ser1NoADC;
ber1NoADC_48_100(ber1NoADC_48_100 > ser_th) = 1;
effNoDiff = pilotUsedChannelEst/Number_of_subcarriers;
throughputNoADC_48_100_2 = (effNoDiff) *(1-ber1NoADC_48_100) * bandwidth;


U = 64;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADCSer/" +strNoADC+ ".mat");
ber1NoADC_64_100 = ber1bitU128_.ser1NoADC;
ber1NoADC_64_100(ber1NoADC_64_100 > ser_th) = 1;
effNoDiff = pilotUsedChannelEst/Number_of_subcarriers;
throughputNoADC_64_100_2 = (effNoDiff) *(1-ber1NoADC_64_100) * bandwidth;



trmsValue = 150;
U = 48;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADCSer/" +strNoADC+ ".mat");
ber1NoADC_48_150 = ber1bitU128_.ser1NoADC;
ber1NoADC_48_150(ber1NoADC_48_150 > ser_th) = 1;
effNoDiff = pilotUsedChannelEst/Number_of_subcarriers;
throughputNoADC_48_150_2 = (effNoDiff) *(1-ber1NoADC_48_150) * bandwidth;



U = 64;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADCSer/" +strNoADC+ ".mat");
ber1NoADC_64_150 = ber1bitU128_.ser1NoADC;
ber1NoADC_64_150(ber1NoADC_64_150 > ser_th) = 1;
effNoDiff = pilotUsedChannelEst/Number_of_subcarriers;
throughputNoADC_64_150_2 = (effNoDiff) *(1-ber1NoADC_64_150) * bandwidth;




%%
end_= 4;
max_x = 20;
max_y = 4;

%ticky = logspace(-6, 0, 7);
%tickx  = linspace(0,max_x,11);
mark_thick = 10;
line_thick =2;
% mark_thick = 5;
figure;
%figure('WindowState','maximized');

%semilogy(SNRdB(1:end_), ber1bitU16_(1:end_),'k-+','LineWidth',line_thick,'MarkerSize',mark_thick)
plot(SNRdB(1:end_), throughputNoDiff_48_50(1:end_),'b--s','LineWidth',line_thick,'MarkerSize',mark_thick)


hold on
grid on; 
plot(SNRdB(1:end_), throughputNoDiff_64_50(1:end_),'b-.o','LineWidth',line_thick,'MarkerSize',mark_thick)
%plot(SNRdB(1:end_), throughputNoDiff_48_100(1:end_),'b-+','LineWidth',line_thick,'MarkerSize',mark_thick)
%plot(SNRdB(1:end_), throughputNoDiff_64_100(1:end_),'b-p','LineWidth',line_thick,'MarkerSize',mark_thick)
plot(SNRdB(1:end_), throughputNoDiff_48_150(1:end_),'b:>','LineWidth',line_thick,'MarkerSize',mark_thick)
plot(SNRdB(1:end_), throughputNoDiff_64_150(1:end_),'b-h','LineWidth',line_thick,'MarkerSize',mark_thick)

plot(SNRdB(1:end_), throughputNoADC_48_50(1:end_),'r--s','LineWidth',line_thick,'MarkerSize',mark_thick)

plot(SNRdB(1:end_), throughputNoADC_64_50(1:end_),'r-.o','LineWidth',line_thick,'MarkerSize',mark_thick)
%plot(SNRdB(1:end_), throughputNoADC_48_100(1:end_),'r-+','LineWidth',line_thick,'MarkerSize',mark_thick)
%plot(SNRdB(1:end_), throughputNoADC_64_100(1:end_),'r-p','LineWidth',line_thick,'MarkerSize',mark_thick)
plot(SNRdB(1:end_), throughputNoADC_48_150(1:end_),'r:>','LineWidth',line_thick,'MarkerSize',mark_thick)
plot(SNRdB(1:end_), throughputNoADC_64_150(1:end_),'r-h','LineWidth',line_thick,'MarkerSize',mark_thick)

plot(SNRdB(1:end_), throughputNoADC_48_50_2(1:end_),'g--s','LineWidth',line_thick,'MarkerSize',mark_thick)

plot(SNRdB(1:end_), throughputNoADC_64_50_2(1:end_),'g-.o','LineWidth',line_thick,'MarkerSize',mark_thick)
%plot(SNRdB(1:end_), throughputNoADC_48_100_2(1:end_),'g-+','LineWidth',line_thick,'MarkerSize',mark_thick)
%plot(SNRdB(1:end_), throughputNoADC_64_100_2(1:end_),'g-p','LineWidth',line_thick,'MarkerSize',mark_thick)
plot(SNRdB(1:end_), throughputNoADC_48_150_2(1:end_),'g:>','LineWidth',line_thick,'MarkerSize',mark_thick)
plot(SNRdB(1:end_), throughputNoADC_64_150_2(1:end_),'g-h','LineWidth',line_thick,'MarkerSize',mark_thick)

%{
semilogy(SNRdB, ber1_96,'k-p','LineWidth',line_thick,'MarkerSize',mark_thick)

semilogy(SNRdB, ber1ADC_16,'r-s','LineWidth',line_thick,'MarkerSize',mark_thick)


semilogy(SNRdB, ber1ADC_48,'r-o','LineWidth',line_thick,'MarkerSize',mark_thick)

semilogy(SNRdB, ber1ADC_96,'r-p','LineWidth',line_thick,'MarkerSize',mark_thick)

semilogy(SNRdB, UnAttber1ADC_16,'b-s','LineWidth',line_thick,'MarkerSize',mark_thick)


semilogy(SNRdB, UnAttber1ADC_48,'b-o','LineWidth',line_thick,'MarkerSize',mark_thick)

semilogy(SNRdB, UnAttber1ADC_96,'b-p','LineWidth',line_thick,'MarkerSize',mark_thick)
%}
%Plot
hold off;

set(gca,'FontSize',25)

xlabel('SNR (dB)' , 'FontSize', 30, 'Interpreter', 'latex')
ylabel('Spectral Efficiency (b/s/Hz)', 'FontSize', 30, 'Interpreter', 'latex')
legend({'Diff: $L = 11$ , $U = 48$';'Diff: $L = 11$ , $U = 64$';'Diff: $L = 31$ , $U = 48$';'Diff: $L = 31$ , $U = 64$';...
'$\xi = 0.125$, $L = 11$ , $U = 48$';'$\xi = 0.125$, $L = 11$ , $U = 64$';'$\xi = 0.125$, $L = 31$ , $U = 48$';'$\xi = 0.125$, $L = 31$ , $U = 64$';    
'$\xi = 0.5$, $L = 11$ , $U = 48$';'$\xi = 0.5$, $L = 11$ , $U = 64$';'$\xi = 0.5$, $L = 31$ , $U = 48$';'$\xi = 0.5$, $L = 31$ , $U = 64$';    
},'Location', 'northwest', 'FontSize', 24, 'Interpreter', 'latex')
%yticks(ticky)
%xticks(tickx)
%xlim([0, max_x])
%ylim([  0, max_y])

%legend({'OSTBC'; 'NN: 100 symbols'; 'NN: 200 symbols'; 'NN: 2000 symbols'; 'NN: 20000 symbols';'NN: 20000, $\rho = 0.7$'}, 'Location', 'southwest', 'FontSize', 35, 'Interpreter', 'latex')
%title('Rate: 2 bit/s/Hz  ', 'FontSize', 12, 'Interpreter', 'latex')


