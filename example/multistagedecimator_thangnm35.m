Astop = 80;
%overall decimation factor of M
M = 16;

%sampling frequency
Fs = 6000;
TW = 0.03*Fs/2;
cMinCoeffs = designMultistageDecimator(M,Fs,TW,Astop,'MinTotalCoeffs',true);
generateFilteringCode(cMinCoeffs,'myDecimator');


%% Example MATLAB code for Multi-Stage Decimation and Bit-Width Optimization
% 
% % Sampling rate and desired decimation factor
% Fs = 10e6;   % Original sampling rate: 10 MHz
% decimation_factors = [5, 5, 2];  % Decimation factors for each stage
% 
% % Stage 1: Low-pass filter design (Decimate by 5)
% Fc1 = 2e6;  % Cutoff frequency for Stage 1 filter (2 MHz)
% order1 = 100;  % Filter order
% b1 = fir1(order1, Fc1/(Fs/2));  % FIR filter design
% 
% % Stage 2: Low-pass filter design (Decimate by 5)
% Fc2 = 0.8e6;  % Cutoff frequency for Stage 2 filter (0.8 MHz)
% order2 = 100;  % Filter order
% b2 = fir1(order2, Fc2/(Fs/2));  % FIR filter design
% 
% % Stage 3: Low-pass filter design (Decimate by 2)
% Fc3 = 0.15e6;  % Cutoff frequency for Stage 3 filter (0.15 MHz)
% order3 = 100;  % Filter order
% b3 = fir1(order3, Fc3/(Fs/2));  % FIR filter design
% 
% % Simulate the decimation process
% % Create a sample signal (e.g., a sine wave at 1 MHz)
% t = 0:1/Fs:1-1/Fs;
% x = sin(2*pi*1e6*t);  % 1 MHz sine wave signal
% 
% % Stage 1: Apply filter and decimate
% x1 = filter(b1, 1, x);  % Apply Stage 1 filter
% x1_down = downsample(x1, 5);  % Decimate by 5
% 
% % Stage 2: Apply filter and decimate
% x2 = filter(b2, 1, x1_down);  % Apply Stage 2 filter
% x2_down = downsample(x2, 5);  % Decimate by 5
% 
% % Stage 3: Apply filter and decimate
% x3 = filter(b3, 1, x2_down);  % Apply Stage 3 filter
% x3_down = downsample(x3, 2);  % Decimate by 2
% 
% % Plot the original and decimated signals
% figure;
% subplot(3,1,1);
% plot(t, x); title('Original Signal');
% subplot(3,1,2);
% plot(x1_down); title('After Stage 1 Decimation');
% subplot(3,1,3);
% plot(x3_down); title('After All Stages Decimation');