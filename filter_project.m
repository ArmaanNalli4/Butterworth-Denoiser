%1) Audio Frequency Analysis

%signal is an array containing the audio data
%Fs is the sampling frequency in Hz (Samples per second)
[signal, Fs] = audioread('noisyaudio.wav');

disp (Fs); %display the sampling frequency

%using the fast fourier transform function to calculate the DFT
N = length (signal);
x = fft(signal);
f = (0:N-1)*(Fs/N); %frequency axis in Hz
mag_dB = 20*log10(abs(x));

% Plot only positive frequencies
half = 1:floor(N/2);

% Plot magnitude in dB
figure;
plot(f(half), mag_dB(half));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Magnitude Spectrum of Noisy Audio (in dB)');
grid on;

hold on;
xline(100, 'g--', '100 Hz');
xline(2000, 'g--', '2 kHz');
xline(2500, 'r--', '2.5 kHz');
xline(5500, 'r--', '5.5 kHz');
legend('Magnitude (dB)', 'Speech Band', 'Noise Band');

%2) Filter Design 
T = 1/Fs; %sampling period

fp = 2000; %Passband edge Frequency (Hz)
fs = 2500; %Stopband edge Frequency (Hz)

%digital passband edge (rad/sample)
w_p = (2* pi * fp)/Fs;
w_s = (2* pi * fs)/Fs;

%converting to analog freq using impulse variance 
W_p = w_p / T;
W_s = w_s / T;

%passband and stopband attenuation
d_p = 10^(-1/20); %passband attenuation <= 1 dB
d_s = 10^(-50/20); %stopband attenuation => 50 dB

%intermediate values for Butterworth eq
kp = (1/d_p^2)-1;
ks = (1/d_s^2)-1;

% Calculate the minimum filter order
N = ceil(0.5 * log10(ks/kp) / log10(W_s/W_p));

% Display filter order
fprintf('Filter order N = %d\n', N);

% Calculate Butterworth cutoff frequency (analog)
W_c = W_p / (kp^(1/(2*N)));
fprintf('Analog cutoff frequency Ω_c = %.4f rad/s\n', W_c);

% Equivalent digital cutoff frequency
w_c = W_c * T;
fprintf('Digital cutoff frequency ω_c = %.4f rad/sample\n', w_c);

% Frequency range for plotting analog gain
Omega = linspace(0, 3*W_c, 1000);

% Butterworth gain in dB
H_dB = -10 * log10(1 + (Omega./W_c).^(2*N));

% Plot gain
figure;
plot(Omega, H_dB, 'LineWidth', 1.5);
xlabel('\Omega (rad/s)');
ylabel('|H_a(j\Omega)| (dB)');
title('Analog Butterworth Filter Logarithmic Gain');
grid on;
ylim([-100 5]);
xline(W_c, 'r--', 'Cutoff \Omega_c');

%Filter Implementation 

Wn = w_c / pi; %nomrmalizing w_c by pi

%Digital Butterworth Low-Pass Filter
[b,a] = butter(N,Wn);

filter_signal = filter(b,a ,signal);

Nf = length(filter_signal);
Xf = fft(filter_signal);
f_axis = (0:Nf-1)*(Fs/Nf);
mag_dB_filtered = 20*log10(abs(Xf));

% Plot magnitude spectrum of filtered audio
half = 1:floor(Nf/2);
figure;
plot(f_axis(half), mag_dB_filtered(half));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Magnitude Spectrum of Filtered Audio (in dB)');
grid on;
hold on;
xline(2000, 'g--', '2 kHz');
xline(2500, 'r--', '2.5 kHz');
xline(5500, 'r--', '5.5 kHz');
legend('Filtered Magnitude (dB)', 'Passband', 'Stopband');

% Compare unfiltered and filtered spectra

figure;
plot(f(half), mag_dB(half), 'b', 'DisplayName', 'Unfiltered');
hold on;
plot(f_axis(half), mag_dB_filtered(half), 'r', 'DisplayName', 'Filtered');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Comparison of Unfiltered vs Filtered Audio Spectrum');
grid on;
legend;

%orignal audio
disp('Playing original audio...');
sound(signal, Fs);
pause(length(signal)/Fs + 1);

%filtered audio
disp('Playing filtered audio...');
sound(filter_signal, Fs);

audiowrite('filteredaudio.wav', filter_signal, Fs);
disp('Filtered audio written to filteredaudio.wav');