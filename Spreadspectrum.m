clc; clear; close all;

% ==========================
% 1. Constants & Parameters
% ==========================
fs = 96000;       % Sampling rate (Hz)
T = 10;           % Duration (s)
order = 7;       % Order for MLS
chip_rate = 50000;  % Chip rate (Hz)
N = fs * T;       % Total number of samples
L = 2^order - 1;  % MLS sequence length
c = 343;          % Speed of sound (m/s)

% Propagation distances
d1 = 100.0;  % Direct path (m)
d2 = 100.3;  % Reflected path (m)

% Define octave band frequencies for visualization
octave_frequencies = [16, 31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000];  

% Time Vector
t = (0:N-1) / fs;

% ==========================
% 2. Generate MLS Sequence (Using Custom Function) ✅
% ==========================
mls_sequence = generate_mls(order);  % Call custom MLS function

% Plot the Generated MLS Sequence
figure;
plot(mls_sequence, 'b');
xlabel('Sample Index');
ylabel('Amplitude');
title('Generated MLS Sequence');
grid on;
ylim([-1.2, 1.2]);

% ==========================
% 3. Expand MLS for Chip Rate ✅
% ==========================
samples_per_chip = round(fs / chip_rate);  % Number of samples per chip
expanded_mls = repelem(mls_sequence, samples_per_chip);  % Stretch MLS sequence

% Ensure the expanded signal matches the desired length
MLS_repeated = repmat(expanded_mls, 1, ceil(N / length(expanded_mls)));
MLS_repeated = MLS_repeated(1:N); % Trim to exact length

% ==========================
% 4. Generate Carrier Wave ✅
% ==========================
carrier = sin(2 * pi * 1000 * t);  % Carrier frequency of 1000 Hz

% ==========================
% 5. Modulate MLS with Carrier ✅
% ==========================
modulated_signal = MLS_repeated .* carrier;  % BPSK Modulation

% ==========================
% 6. Compute FFT and Phase Shifts ✅
% ==========================
S = fft(modulated_signal, N);  % Compute FFT of original signal
frequencies = (0:N/2-1) * (fs / N);  % Frequency axis (only positive)

% Compute exact phase shift for each frequency component
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

% ==========================
% 7. Compute Transfer Function ✅
% ==========================
H_f = X(1:N/2) ./ (S(1:N/2) + eps);  % H(f) = X(f) / S(f)
H_mag_dB = 20 * log10(abs(H_f) + eps);  % Magnitude in dB
H_phase = unwrap(angle(H_f)) * (180/pi);  % Phase response in degrees

% ==========================
% 8. Plot Transfer Function ✅
% ==========================
figure;

subplot(2,1,1);
semilogx(frequencies, H_mag_dB);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Transfer Function Magnitude (dB)');
xlim([50,9000])
grid on; hold on;

% Add vertical lines at octave band frequencies
for f = octave_frequencies
    xline(f, '--r', 'LineWidth', 1.2);
end
hold off;

subplot(2,1,2);
semilogx(frequencies, H_phase);
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
title('Transfer Function Phase Response');
grid on; hold on;

% Add vertical lines at octave band frequencies
for f = octave_frequencies
    xline(f, '--r', 'LineWidth', 1.2);
end
hold off;

% ==========================
% 9. Define Custom MLS Function ✅
% ==========================
function seq = generate_mls(m)
    N = 2^m - 1;
    reg = ones(1, m);  % Initial state (vector of 1s)
    seq = zeros(1, N); % Preallocate output sequence

    for i = 1:N
        seq(i) = reg(end);  % Output bit
        feedback = xor(reg(1), reg(end));  % XOR feedback
        reg = [feedback, reg(1:end-1)];  % Shift register
    end

    seq = 2 * seq - 1;  % Convert {0,1} to {-1,1}
end
