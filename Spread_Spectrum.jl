using Random, FFTW, Plots, StatsBase, DSP
include("MLS.jl")
plotly()

# ==========================
# 1. Constants & Parameters
# ==========================
fs = 48000       # Sampling rate (Hz)
T = 10           # Duration (s)
order = 6       # MLS Order
chip_rate = 200 # Chip rate (Hz)
N = Int(fs * T)  # Total number of samples
f = 1000         # Carrier frequency (Hz)
d1, d2 = 100, 100.3
c = 343

# Time Delays
time_delay_1 = d1 / c  
time_delay_2 = d2 / c  

# Sample Delays
sample_delay_1 = round(Int, time_delay_1 * fs)  
sample_delay_2 = round(Int, time_delay_2 * fs)

# Time Vector
time = range(0, stop=T, length=N)

# -------------------------------
#  Generate MLS Sequence

mls_signal = generate_mls(order)

# 
# ----------------------------------------------
samples_per_chip = round(Int, fs / chip_rate)  
expanded_mls = repeat(mls_signal, inner=samples_per_chip)  

# Ensure correct length
MLS_repeated = repeat(expanded_mls, ceil(Int, N / length(expanded_mls)))[1:N]

# ==========================
#
#------------------MODULATIOIN---------------------
carrier = sin.(2π * f .* time)
spreadsignal = MLS_repeated .* carrier

# Introduce Delays
signal_d1 = spreadsignal  
#signal_d2 = vcat(zeros(sample_delay_1), signal_d1)[1:length(signal_d1)]
signal_d2 = circshift(spreadsignal, sample_delay_2)

# Sum the delayed signals
sumsignal = signal_d1 .+ signal_d2

# ==========================
#FFT
# ==========================
X= fft(sumsignal)
S = fft(spreadsignal)  # Use original signal as reference

# Compute Transfer Function
H_mag = X ./ S
H_mag_dB = 20 * log10.(abs.(H_mag[1:N÷2]))

# Correct Frequency Axis
frequencies = (0:N÷2-1) .* (fs / N)  # Now correctly scaled


# ==========================
# 6. Compute Autocorrelation
# ==========================
AC_nonfiltered = autocor(spreadsignal .- mean(spreadsignal), 0:N÷2)
AC_nonfiltered /= maximum(abs.(AC_nonfiltered))

# ==========================
# 7. Plot Results
# ==========================
p1 = plot(time, sumsignal, xlabel="Time [s]", ylabel="Amplitude", title="Spread Spectrum Signal", label="Modulated MLS")
p2 = plot(mls_signal, xlabel="Samples", ylabel="Amplitude", title="MLS Sequence", label="MLS")
p3 = plot(frequencies, H_mag_dB, xlabel="Frequency [Hz]", ylabel="Magnitude (dB)", title="Transfer Function", xscale=:log10, label="H(f)")
p4 = plot(AC_nonfiltered, xlabel="Lag", ylabel="Autocorrelation", title="Autocorrelation", label="Autocorrelation")

plot(p1, p2, p3, p4, layout=(4,1), size=(1200,1200))
