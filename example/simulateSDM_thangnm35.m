%simulate CIFF for filter design
%OSR = 16
%
% H = synthesizeNTF(5,32,1);
% [a,b,c,d]=realizeNTF(H,'CIFF')

%{
ntf = zpk([1 1], [1 1]/3,1,1)
n_imp = 10;
y_desired = impL1(ntf,n_imp)';
Ac = [0 0; 1 0];
Bc = [-1 0 0; 0 -1 0];
Cc = [0 1];
Dc = [0 0 -1];
td = 0;
sys_c = ss(Ac,Bc,Cc,Dc);
set(sys_c,'InputDelay',td*[1 1 1]);
sys_d = c2d(sys_c,1)
yy = squeeze(impulse(sys_d,n_imp))';
a = y_desired / yy

%}
% 
clear all
close all
%% Design parameters
%filter order
order = 4;
%oversampling rate
osr = 16;
%number of levels
nlev = 16; 
f0 = 0;
% %The maximum out-of-band gain of the NTF
Hinf = 2.5; % 3.162 
% % 2.5 for CT in cadence
% % Parameters for the continuous-time implementation
% tdac = [0 1];	% DAC timing. [0 1] means zero-delay non-return-to-zero
% form = 'FF';
% 
% % Derived parameters
% M = nlev-1;
% clc
% fprintf(1,'\t\t\t%dth-Order Continuous-Time Lowpass Example\n\n',order);
% fig_pos = { [  10 595  600 215],
%             [ 640 595  600 215],
%             [  10 345  300 205],
%             };
% 
% % NTF synthesis and realization
% fprintf(1,'Doing NTF synthesis... ');
% design_step = 1;
% 
% Fs = 3.2e9;
% BW = 100e6;
ntf0 = synthesizeNTF(order,osr,0,Hinf,f0);		% Optimized zero placement
% 
% % 0 is no zero optimization
% % 1 is zero optimization
% 
% fprintf(1,'Done.\n');
% figure(design_step); clf;
% %set(design_step,'position',fig_pos{design_step});
% ntf_axes = DocumentNTF(ntf0,osr,f0); %Plot the NTF's poles and zeros as well as its frequency-response
% drawnow;
% 
% % Time-domain simulations
% fprintf(1,'Doing time-domain simulations... ');
% design_step = design_step+1;
% figure(design_step); clf;
% %set(design_step,'position',fig_pos{design_step});
% % Example spectrum
% subplot('position', [0.05,0.1,0.6 0.8]);
% PlotExampleSpectrum(ntf0,M,osr,f0);
% title('Example Spectrum');
% % SQNR plot
% subplot('position',[.74 .18 .25 .65]);
% if nlev==2
%     [snr_pred,amp_pred] = predictSNR(ntf0,osr);
%     plot(amp_pred,snr_pred,'-');
%     hold on;
% end
% [snr,amp] = simulateSNR(ntf0,osr,[],f0,nlev);
% fprintf(1,'Done.\n');
% plot(amp,snr,'og');
% figureMagic([-100 0], 10, 2, [0 100], 10, 2, [7 3], 'Discrete-Time Simulation');
% xlabel('Input Level (dBFS)');
% ylabel('SQNR (dB)');
% 
% %Find the snr peak by fitting
% [peak_snr,peak_amp] = peakSNR(snr,amp);
% msg = sprintf('peak SQNR = %4.1fdB  \n@ amp=%4.1fdB  ',peak_snr,peak_amp);
% text(peak_amp-10,peak_snr,msg,'hor','right', 'vertical','middle');
% msg = sprintf('OSR=%d ',osr);
% text(0,5,msg,'hor','right');
% title('SQNR Plot');
% drawnow;
% 
% % Continuous-Time Mapping
% fprintf(1,'Mapping  to continuous-time... ');
% design_step = design_step+1;
% %Realize an NTF with a continuous-time loop filter.
% [ABCDc,tdac2] = realizeNTF_ct( ntf0, form, tdac);
% %%% adjust D = b4 to see what happen in STF
% ABCDc(4,5) = 0.00;
% %ABCDc(2,5) = -0.05;
% [Ac Bc Cc Dc] = partitionABCD(ABCDc);
% fprintf(1,'Done.\n');
% 
% 

%% thangnm35
% %book example
% h_ntf = synthesizeNTF(5,32,1);
% N = 8192; f = 85; fB = ceil(N/(2*osr));
% u = 0.5*sin(2*pi*fB/N*[0:N-1]);
% v = simulateDSM(u,h_ntf);
% t = 0:85;
% stairs(t, u(t+1),'g');
% hold on;
% stairs(t,v(t+1),'b');
% axis([0 85 -1.2 1.2]);
% ylabel('u, v');
% 
% figure();
% spec=fft(v.*ds_hann(N))/(N/4);
% plot(linspace(0,0.5,N/2+1), ...
% dbv(spec(1:N/2+1)));
% axis([0 0.5 -120 0]);
% grid on;
% ylabel('dBFS/NBW')
% snr=calculateSNR(spec(1:fB),f);
% s=sprintf('SNR = %4.1fdB\n',snr)
% text(0.25,-90,s);
% s=sprintf('NBW=%7.5f',1.5/N);
% text(0.25, -110, s);

% % figure(); clf;
% Fb=Fs/(2*OSR);  
f_data = 1e3;
mul_factor = 1.1;
f_data_max = f_data*mul_factor; 
N = 2^(ceil(log2(f_data_max*(2*osr))))

user_data = 0.5*sin(2*pi*f_data/N*(0:N-1));
[v,xn,xmax,y] = simulateDSM(user_data, ntf0, nlev, 0);

t = 0:1000;
stairs(t, user_data(t+1),'g');
hold on;
stairs(t, xn(t+1),'b');
axis([0 1000 -1.2 1.2]);
ylabel('user_data, sdm_data');