clc
clear
close all

Nt = 2;
Nr = 1;
SNRdB = linspace(0,15,4);%linspace(0,40,21);
SNRdB = linspace(-5,15,5);%linspace(0,40,21);

SNRList = 10.^(SNRdB/10);


ModulationType = '8-PSK';
M = 8; 
bitsPerSymbol = log2(M);  

fd = 5;
trmsValue = 50;
U = 48;
alowhigh = 2;
formatSpec = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d';
str = sprintf(formatSpec,M,Nt,Nr,U,fd,trmsValue);
ber1bitU128_ = load("Data/ADCRed/" +str+ ".mat");
ber1bitRed_48_50 = ber1bitU128_.ber1bitRed;



U = 64;
alowhigh = 2;
formatSpec = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d';
str = sprintf(formatSpec,M,Nt,Nr,U,fd,trmsValue);
ber1bitU128_ = load("Data/ADCRed/" +str+ ".mat");
ber1bitRed_64_50 = ber1bitU128_.ber1bitRed;


trmsValue = 100;
U = 48;
alowhigh = 2;
formatSpec = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d';
str = sprintf(formatSpec,M,Nt,Nr,U,fd,trmsValue);
ber1bitU128_ = load("Data/ADCRed/" +str+ ".mat");
ber1bitRed_48_100 = ber1bitU128_.ber1bitRed;



U = 64;
alowhigh = 2;
formatSpec = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d';
str = sprintf(formatSpec,M,Nt,Nr,U,fd,trmsValue);
ber1bitU128_ = load("Data/ADCRed/" +str+ ".mat");
ber1bitRed_64_100 = ber1bitU128_.ber1bitRed;



trmsValue = 100;
U = 48;
alowhigh = 2;
formatSpec = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d';
str = sprintf(formatSpec,M,Nt,Nr,U,fd,trmsValue);
ber1bitU128_ = load("Data/ADCRed/" +str+ ".mat");
ber1bitRed_48_150 = ber1bitU128_.ber1bitRed;



U = 64;
alowhigh = 2;
formatSpec = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d';
str = sprintf(formatSpec,M,Nt,Nr,U,fd,trmsValue);
ber1bitU128_ = load("Data/ADCRed/" +str+ ".mat");
ber1bitRed_64_150 = ber1bitU128_.ber1bitRed;



pilotFrac =4;
NtNoADC =1;
fd = 5;
trmsValue = 50;
U = 48;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADC/" +strNoADC+ ".mat");
ber1NoADC_48_50 = ber1bitU128_.ber1NoADC;



U = 64;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADC/" +strNoADC+ ".mat");
ber1NoADC_64_50 = ber1bitU128_.ber1NoADC;


trmsValue = 100;
U = 48;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADC/" +strNoADC+ ".mat");
ber1NoADC_48_100 = ber1bitU128_.ber1NoADC;



U = 64;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADC/" +strNoADC+ ".mat");
ber1NoADC_64_100 = ber1bitU128_.ber1NoADC;



trmsValue = 100;
U = 48;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADC/" +strNoADC+ ".mat");
ber1NoADC_48_150 = ber1bitU128_.ber1NoADC;



U = 64;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADC/" +strNoADC+ ".mat");
ber1NoADC_64_150 = ber1bitU128_.ber1NoADC;



pilotFrac =2;
NtNoADC =1;
fd = 5;
trmsValue = 50;
U = 48;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADC/" +strNoADC+ ".mat");
ber1NoADC_48_50_2 = ber1bitU128_.ber1NoADC;



U = 64;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADC/" +strNoADC+ ".mat");
ber1NoADC_64_50_2 = ber1bitU128_.ber1NoADC;


trmsValue = 100;
U = 48;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADC/" +strNoADC+ ".mat");
ber1NoADC_48_100_2 = ber1bitU128_.ber1NoADC;



U = 64;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADC/" +strNoADC+ ".mat");
ber1NoADC_64_100_2 = ber1bitU128_.ber1NoADC;



trmsValue = 100;
U = 48;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADC/" +strNoADC+ ".mat");
ber1NoADC_48_150_2 = ber1bitU128_.ber1NoADC;



U = 64;
alowhigh = 2;
formatSpecNoADC = '_M_%d_Nt_%d_Nr_%d_U_%d_fd_%d_trms_%d_pilot_%d';
strNoADC = sprintf(formatSpecNoADC,M,NtNoADC,Nr,U,fd,trmsValue,pilotFrac);
ber1bitU128_ = load("Data/NoADC/" +strNoADC+ ".mat");
ber1NoADC_64_150_2 = ber1bitU128_.ber1NoADC;




%%
end_= 4;
max_x = 20;
max_y = 10^-5;

%ticky = logspace(-6, 0, 7);
%tickx  = linspace(0,max_x,11);
mark_thick = 10;
line_thick = 3;
% mark_thick = 5;
figure;

%figure('WindowState','maximized');

%semilogy(SNRdB(1:end_), ber1bitU16_(1:end_),'k-+','LineWidth',line_thick,'MarkerSize',mark_thick)
semilogy(SNRdB(1:end_), ber1bitRed_48_50(1:end_),'b--s','LineWidth',line_thick,'MarkerSize',mark_thick)


hold on
grid on; 
semilogy(SNRdB(1:end_), ber1bitRed_64_50(1:end_),'b--o','LineWidth',line_thick,'MarkerSize',mark_thick)
%semilogy(SNRdB(1:end_), ber1bitRed_48_100(1:end_),'b-+','LineWidth',line_thick,'MarkerSize',mark_thick)
%semilogy(SNRdB(1:end_), ber1bitRed_64_100(1:end_),'b-p','LineWidth',line_thick,'MarkerSize',mark_thick)
semilogy(SNRdB(1:end_), ber1bitRed_48_150(1:end_),'b->','LineWidth',line_thick,'MarkerSize',mark_thick)
semilogy(SNRdB(1:end_), ber1bitRed_64_150(1:end_),'b-h','LineWidth',line_thick,'MarkerSize',mark_thick)

semilogy(SNRdB(1:end_), ber1NoADC_48_50(1:end_),'r--s','LineWidth',line_thick,'MarkerSize',mark_thick)

semilogy(SNRdB(1:end_), ber1NoADC_64_50(1:end_),'r--o','LineWidth',line_thick,'MarkerSize',mark_thick)
%semilogy(SNRdB(1:end_), ber1NoADC_48_100(1:end_),'r-+','LineWidth',line_thick,'MarkerSize',mark_thick)
%semilogy(SNRdB(1:end_), ber1NoADC_64_100(1:end_),'r-p','LineWidth',line_thick,'MarkerSize',mark_thick)
semilogy(SNRdB(1:end_), ber1NoADC_48_150(1:end_),'r->','LineWidth',line_thick,'MarkerSize',mark_thick)
semilogy(SNRdB(1:end_), ber1NoADC_64_150(1:end_),'r-h','LineWidth',line_thick,'MarkerSize',mark_thick)

semilogy(SNRdB(1:end_), ber1NoADC_48_50_2(1:end_),'g--s','LineWidth',line_thick,'MarkerSize',mark_thick)

semilogy(SNRdB(1:end_), ber1NoADC_64_50_2(1:end_),'g--o','LineWidth',line_thick,'MarkerSize',mark_thick)
%semilogy(SNRdB(1:end_), ber1NoADC_48_100_2(1:end_),'r-+','LineWidth',line_thick,'MarkerSize',mark_thick)
%semilogy(SNRdB(1:end_), ber1NoADC_64_100_2(1:end_),'r-p','LineWidth',line_thick,'MarkerSize',mark_thick)
semilogy(SNRdB(1:end_), ber1NoADC_48_150_2(1:end_),'g->','LineWidth',line_thick,'MarkerSize',mark_thick)
semilogy(SNRdB(1:end_), ber1NoADC_64_150_2(1:end_),'g-h','LineWidth',line_thick,'MarkerSize',mark_thick)

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
ylabel('BER', 'FontSize', 30, 'Interpreter', 'latex')
legend({'Diff: $L = 11$ , $U = 48$';'Diff: $L = 11$ , $U = 64$';'Diff: $L = 31$ , $U = 48$';'Diff: $L = 31$ , $U = 64$';...
'$\xi = 0.5$, $L = 11$ , $U = 48$';'$\xi = 0.5$, $L = 11$ , $U = 64$';'$\xi = 0.5$, $L = 31$ , $U = 48$';'$\xi = 0.5$, $L = 31$ , $U = 64$';    
'$\xi = 0.125$, $L = 11$ , $U = 48$';'$\xi = 0.125$, $L = 11$ , $U = 64$';'$\xi = 0.125$, $L = 31$ , $U = 48$';'$\xi = 0.125$, $L = 31$ , $U = 64$';    
},'Location', 'southwest', 'FontSize', 24, 'Interpreter', 'latex')
%yticks(ticky)
%xticks(tickx)
%xlim([0, max_x])
%ylim([  0, max_y])
ylim([max_y 1 ])

axes('position',[.45 .175 .25 .25])
box on % put box around new pair of axes
indexOfInterest = (SNRdB <= 10) & (SNRdB >= 5); % range of t near perturbation
semilogy(SNRdB(indexOfInterest), ber1bitRed_48_50(indexOfInterest),'b--s','LineWidth',line_thick,'MarkerSize',mark_thick)
hold on
semilogy(SNRdB(indexOfInterest), ber1bitRed_48_150(indexOfInterest),'b->','LineWidth',line_thick,'MarkerSize',mark_thick)
hold off
axis tight


axes('position',[.45 .175 .25 .25])
box on % put box around new pair of axes
indexOfInterest = (SNRdB <= 10) & (SNRdB >= 5); % range of t near perturbation
semilogy(SNRdB(indexOfInterest), ber1bitRed_64_50(indexOfInterest),'b--o','LineWidth',line_thick,'MarkerSize',mark_thick)
hold on
semilogy(SNRdB(indexOfInterest), ber1bitRed_64_150(indexOfInterest),'b-h','LineWidth',line_thick,'MarkerSize',mark_thick)
hold off
axis tight

