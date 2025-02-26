% Parameters
order = 5
L = 2^order-1; % Order of the MLS sequence
fs = 48000; % Sampling frequency (Hz)
duration = 10; % Duration of the signal in seconds

% Generate MLS sequence
mls_sequence = mls(order); % Generates an MLS sequence of length 2^order - 1

% Calculate the required length to match the white noise sequence
required_length = fs * duration;

% Repeat the MLS sequence to match the required length
repeated_mls = repmat(mls_sequence, ceil(required_length / length(mls_sequence)), 1);
repeated_mls = repeated_mls(1:required_length); % Truncate to match exact length

% Generate white noise sequence
white_noise = randn(required_length, 1); % Gaussian white noise

% Normalize both sequences to have the same RMS level
repeated_mls = repeated_mls / rms(repeated_mls);
white_noise = white_noise / rms(white_noise);

% Modulate the MLS sequence into the white noise
modulated_signal = white_noise + repeated_mls;

% Compute the FFT of the white noise, MLS sequence, and modulated signal
N = length(white_noise); % Length of the signals
frequencies = (0:N/2-1) * (fs / N); % Frequency axis (first half)

% Compute FFT and take only the first half (positive frequencies)
white_noise_fft = abs(fft(white_noise) / max(abs(fft(white_noise))));
mls_fft = abs(fft(repeated_mls) / max(abs(fft(repeated_mls))));
modulated_fft = abs(fft(modulated_signal) / max(abs(fft(modulated_signal))));

white_noise_fft_half = white_noise_fft(1:N/2);
mls_fft_half = mls_fft(1:N/2);
modulated_fft_half = modulated_fft(1:N/2);

% Time axis for plotting
t = (0:length(repeated_mls)-1) / fs;

% Plot the time-domain signals
figure;
subplot(3,1,1);
plot(t, repeated_mls);
title('Repeated MLS Sequence');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,1,2);
plot(t, white_noise);
title('White Noise Sequence');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,1,3);
plot(t, modulated_signal);
title('Modulated Signal (White Noise + MLS)');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot the FFT of the signals
figure;
subplot(3,1,1);
plot(frequencies, 20*log10(white_noise_fft_half));
title('FFT of White Noise');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
xlim([0 fs/2]);

subplot(3,1,2);
plot(frequencies, 20*log10(mls_fft_half));
title('FFT of MLS Sequence');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
xlim([0 fs/2]);

subplot(3,1,3);
plot(frequencies, 20*log10(modulated_fft_half));
title('FFT of Modulated Signal (White Noise + MLS)');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
xlim([0 fs/2]);
