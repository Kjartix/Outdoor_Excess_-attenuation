clc; clear; close all;

% Define chirp parameters
f_start = 16;   % Start frequency in Hz
f_end = 10000;  % Ending frequency in Hz
T = 10;         % Duration in seconds
Fs = 48000;     % Sampling frequency in Hz
c = 343;        % Speed of sound in m/s
d1 = 100;       % Propagation distance 1 (m)
d2 = 100.3;     % Propagation distance 2 (m)
m = 0:1:7
notch_frequncies = (2*m+1)*c ./(2*(0.3))
% SNR settings

desired_SNR_dB = 7;  % Adjust SNR level (higher = less noise)

% Time vector
t = 0:1/Fs:T;
t = t(1:end-1);  % Ensure consistent length
N = length(t);   % Number of points

% Generate a linear chirp
signal = chirp(t, f_start, T, f_end, 'linear');

% Compute sample delays for each propagation distance
delay_samples_1 = round((d1 / c) * Fs);
delay_samples_2 = round((d2 / c) * Fs);

% Create delayed signals by padding with zeros
signal_d1 = [zeros(1, delay_samples_1), signal(1:end-delay_samples_1)];
signal_d2 = [zeros(1, delay_samples_2), signal(1:end-delay_samples_2)];

% Apply distance-based attenuation
attenuation_1 = 1 / d1;  % Inverse distance attenuation
attenuation_2 = 1 / d2;
signal_d1 = signal_d1 * attenuation_1;
signal_d2 = signal_d2 * attenuation_2;

% Add noise to attenuated signals
signal_d1_power = mean(signal_d1.^2);
signal_d2_power = mean(signal_d2.^2);
noise_d1_power = signal_d1_power / (10^(desired_SNR_dB / 10));
noise_d2_power = signal_d2_power / (10^(desired_SNR_dB / 10));
noise_d1 = randn(1, N) * sqrt(noise_d1_power);
noise_d2 = randn(1, N) * sqrt(noise_d2_power);
signal_d1 = signal_d1 + noise_d1;
signal_d2 = signal_d2 + noise_d2;

% Sum the two signals
x = signal_d1 + signal_d2;

% Plot signals in time domain
figure;
subplot(4,1,1);
plot(t, signal, 'k');
title('Original Chirp Signal');
xlabel('Time (s)'); ylabel('Amplitude');
xlim([0 0.1]);  % Zoom into the first 0.1 seconds
grid on;

subplot(4,1,2);
plot(t, signal_d1, 'b');
title('Signal After 100 m Propagation');
xlabel('Time (s)'); ylabel('Amplitude');
xlim([0 10]);
grid on;

subplot(4,1,3);
plot(t, signal_d2, 'r');
title('Signal After 100.3 m Propagation');
xlabel('Time (s)'); ylabel('Amplitude');
xlim([0 10]);
grid on;

subplot(4,1,4);
plot(t, x, 'm');
title('Summed Signal (100m + 100.3m)');
xlabel('Time (s)'); ylabel('Amplitude');
xlim([0 10]);
grid on;

% Compute FFT of summed signal
X = fft(x, N);

% Compute Transfer Function H(f) = X(f) / S(f)
S = fft(signal, N);  % FFT of the original signal
H_f = X ./ (S + eps);  % Avoid division by zero
H_mag = abs(H_f);  % Magnitude response
H_mag_dB = 20 * log10(H_mag + eps);  % Convert to dB

% Compute frequency vector
frequencies = (0:N/2-1) * (Fs / N);

% Plot Transfer Function Magnitude
figure;
semilogx(frequencies, H_mag_dB(1:N/2), 'b', 'DisplayName', 'Magnitude Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Transfer Function Magnitude using Time-Domain Delays');
xlim([15, 8000]);
grid on; hold on;

for f = notch_frequncies
    xline(f, '--r', 'LineWidth', 1.2);
end
