clc; clear; close all;

% Define chirp parameters
f_start = 16;   % Start frequency in Hz
f_end = 10000;   % Ending frequency in Hz
T = 10;         % Duration in seconds
Fs = 96000;     % Sampling frequency in Hz
c = 343;        % Speed of sound in m/s
d1 = 100;       % Propagation distance 1 (m)
d2 = 100.3;     % Propagation distance 2 (m)
octave_frequencies = [16, 31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000];  % Hz

% Time vector
t = 0:1/Fs:T;
t = t(1:end-1);  % Ensure consistent length
N = length(t);   % Number of points

% Generate a linear chirp
signal = chirp(t, f_start, T, f_end, 'linear');

% Compute FFT of the original signal
S = fft(signal, N);  % Ensure FFT length matches N

% Compute frequency vector (only positive frequencies)
frequencies = (0:N/2-1) * (Fs / N);

% Compute exact phase shift due to delay for each frequency component
time_delay_1 = d1 / c;  
time_delay_2 = d2 / c;  
phase_shift_1 = exp(-1j * 2 * pi * frequencies * time_delay_1);
phase_shift_2 = exp(-1j * 2 * pi * frequencies * time_delay_2);

% Apply phase shift in frequency domain
S_pos = S(1:N/2);  % Take only positive frequencies for symmetry
X1_pos = S_pos .* phase_shift_1;  % Delayed version 1
X2_pos = S_pos .* phase_shift_2;  % Delayed version 2

% Reconstruct full spectrum (using conjugate symmetry)
X1 = [X1_pos, conj(flip(X1_pos(2:end)))];  % Complete with symmetry
X2 = [X2_pos, conj(flip(X2_pos(2:end)))];

% Convert back to time domain (IFFT)
signal_d1 = real(ifft(X1, N));  
signal_d2 = real(ifft(X2, N));

% Sum the two signals
x = signal_d1 + signal_d2;

% Compute FFT of summed signal
X = fft(x, N);  % Ensure FFT size matches N

% Compute Transfer Function H(f) = X(f) / S(f)
H_f = X(1:N/2) ./ (S(1:N/2) + eps);  % Ensure equal sizes
H_mag = abs(H_f);  % Magnitude response
H_mag_dB = 20 * log10(H_mag + eps);  % Convert to dB
H_phase = angle(H_f);  % Phase response in radians

% Create plots
figure;
subplot(4,1,1);
plot(t, x);
xlabel('Time (s)'); ylabel('Amplitude');
title('Summed Signal at 100m and 100.3m');

subplot(4,1,2);
plot(t, signal_d1);
xlabel('Time (s)'); ylabel('Amplitude');
title('Chirp Signal at 100m');

subplot(4,1,3);
plot(t, signal_d2);
xlabel('Time (s)'); ylabel('Amplitude');
title('Chirp Signal at 100.3m');

subplot(4,1,4);
semilogx(frequencies, H_mag_dB);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Transfer Function Magnitude (dB)');
grid on; hold on;

% Add vertical lines at octave band frequencies
for f = octave_frequencies
    xline(f, '--r', 'LineWidth', 1.2);
end
hold off;

% Plot Phase Response
figure;
semilogx(frequencies, rad2deg(H_phase)); % Convert phase to degrees
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
title('Transfer Function Phase Response');
grid on; hold on;

% Add vertical lines at octave band frequencies
for f = octave_frequencies
    xline(f, '--r', 'LineWidth', 1.2);
end
hold off;
