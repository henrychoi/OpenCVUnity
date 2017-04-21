function matched_lfm
musicDir = '~/Music/';
osxDataDir = '~/radar/Qosx/';
iosDataDir = '~/radar/QiOS/';
% Write the waveform and filter for C prog
osxSrcDir = strcat(osxDataDir, 'Qosx/');

TRUE_DISTANCE = 6
AMBIENT_SOUND = 0
CHIRP_STRENGTH = 1

MULTIPATH_DISTANCE = TRUE_DISTANCE + 0.02;
Fs = 44100;
Ts = 1/Fs;
% MacBook Pro fades in on the low end around 130 Hz
f0 = 5000; % lowest frequency of the chirp
b = 10000; % bandwidth
M = 2^13; % # samples to describe the original chirp (controls tau)
N = M * 2; %Total samples I need to collect
c = 335; % speed of sound [m/s]
taup = M * Ts;
Rmax = N * Ts * c;
freqlimit = 0.5 * Fs;

%fprintf('Detection range: %f m\n', 0.5 * Nsample * Ts * c);
fprintf('Range resolution: %f mm\n', 1000*c/b);

mu = b/taup; % compression ratio; @ t=taup, instantaneous F = f0 + b
scat_range = [TRUE_DISTANCE];% MULTIPATH_DISTANCE];
scat_rcs = [1 0.5]; %[1 1.5 2];
nscat = length(scat_range);
%eps = 1.0e-16;

t = Ts * (0:N-1);
t_xmit = t(1:M);

switch(AMBIENT_SOUND)
    case 1
        [s_ambient, Fs_a] = audioread('singing.wav', [1 2*N] + 8E4);
        if Fs_a ~= Fs % resample
            s_ambient = interp1(((1:length(s_ambient)) - 1) / Fs_a, ...
                s_ambient, t, 'linear');%, 'extrap')';
        else
            s_ambient = s_ambient';
        end
    case 2
        s_ambient = sin(10000*2*pi*Ts*(1:N));
    otherwise
        s_ambient = zeros(1,N);
end
RMS_ambient = sqrt(mean(s_ambient .^ 2));
sigma_r = 0.1 * max(RMS_ambient, 1);
A_u = CHIRP_STRENGTH * max(RMS_ambient * sqrt(2), 1);

% time bandwidth product
time_B_product = b * taup;
if(time_B_product < 5 )
    fprintf('************ Time Bandwidth product is TOO SMALL ***************')
    fprintf('\n Change b and or taup')
  return
end
% check to ensure that scatterers are within recieve window
index = find(scat_range > Rmax);
if (index ~= 0)
    error('Error. Receive window is too large; or scatterers fall outside window')
end

% initialize input, output and xmit vectors
x(nscat,1:N) = 0.;
rcv_noise = sigma_r * randn(1,N);
s_r = rcv_noise;
s_r = s_r + s_ambient((1:N));

% Form the transmitted signal
pulse = sin(2*pi * (f0 .* t_xmit + 0.5*mu .* t_xmit.^2) ...
    + 0.1*pi); % Start with non-zero amplitude
xmit = A_u * pulse; % Adjust pulse loudness

silence = zeros(size(pulse));
audiowrite(strcat(musicDir, 'rchirp.wav'), [silence; pulse]', Fs);
csvwrite(strcat(osxSrcDir, 'chirp.h'), pulse);
sound(xmit, Fs);% What does the chirp sound like?

%return;

k_true = round(Fs * scat_range(1) / c);
sdots = dot(xmit, pulse) / (scat_range(1).^2);
sdotn = zeros(1, 20);
for i=1:length(sdotn)
    sdotn(i) = dot(s_r((k_true + i - length(sdotn)/2) + (1:M)), pulse);
end
% figure(4); subplot(212); plot(sdotn)

% Simulate multi-path
for j = 1:nscat
    range = scat_range(j);
    % Discretize the distance into the sample instance
    t_R = range/c; % nominal wave travel time
    x(j, (1:M) + round(t_R * Fs)) = xmit;
    
    % scale receoved signal and accumulate
    s_r = s_r + (scat_rcs(j) / (range)^2) * x(j,:);
end
sdotn = dot(s_r(k_true + (1:M)), pulse);

figure(2);
subplot(211); plot(t_xmit, xmit, t(1:N), s_ambient((1:N)), ':', t(1:N), rcv_noise, '--');
title('truth'); legend('xmit', 'ambient', 'noise');
subplot(212); semilogy(t(1:M), abs(xmit) .^2, t,abs(s_r).^2, ':'); %grid
title('Sound amplitude'); legend('xmit', 'received'); xlabel('[s]');
axis([-inf inf 1E-3, inf])

%% Correlate in time domain
%pulse = pulse(end:-1:1); % reverse it
txrx_xcorr = xcorr(s_r, pulse); % len = 2*N - 1

figure(3);
subplot(411); plot([(t(1:(end-1)) - t(end)) t], txrx_xcorr); title('xcorr');
%rx_pad = zeros(size(pulse_pad) + [1 n]);
% See my blog
% for sliding corr interval for why there are only N-M+1 valid correlations
% But to keep the variables a power of 2, I brazenly throw away the last
% valid correlation, which IS valid.
time_corr = zeros(size(s_r) - [0 M]);
corr = xcorr(s_r(1:M), pulse); % size = 2M-1; 
time_corr(1:M) = corr(M:end); %Only the middle sample is complete;
%the samples to the right have to be saved. M-1 samples on the left are
%invalid for the application because they are incomplete (and remain so).
subplot(413); plot(t(1:length(time_corr)), time_corr);

for i=M:M:(N-2*M) % Feed subsequent received samples
    corr = xcorr(s_r(i+(1:M)), pulse); % dim = 2*M-1
    % Add [1:M-1] samples to the previous corr's [M+1:end] samples to
    % complete the corr, [M] is complete, and [M+1:end] have to be saved
    % for the next round.  So ALL of corr should just be ADDED to txrx_corr
    % because txrx_corr is initialized to 0
    time_corr(i + (-(M-2):M)) = time_corr(i + (-(M-2):M)) + corr;
    % Debug incremental update
    subplot(413); plot(t(1:length(time_corr)), time_corr);
end
% Pick up the last block's contribution (L half: M-1)
if isempty(i), i = M; else i = i + M; end
corr = xcorr(s_r(i+(1:M)), pulse);
time_corr(i + (-(M-2):0)) = time_corr(i + (-(M-2):0)) + corr(1:(M-1));
ha(3) = subplot(413); plot(t(1:length(time_corr)), time_corr);
title('time corr');

% Check my algorithm against Matlab's xcorr
ha(2) = subplot(412);
plot(t(1:length(time_corr)), txrx_xcorr(N-1 + (1:length(time_corr)) ) );
title('time corr should match this (xcorr)');

err = txrx_xcorr(N-1 + (1:length(time_corr))) - time_corr;
if max(abs(err)) > 1E-9
    ha(4) = subplot(414); plot(t(1:length(time_corr)), err);
    error('Overlapped correlation is wrong; e = %f', max(abs(err)))
end

figure(1)
subplot(2,1,1); plot(t(1:M),real(xmit)); grid
xlabel('t [s]'); ylabel('real(xmit)');

freq = linspace(-freqlimit,freqlimit, M);

subplot(2,1,2);
semilogy(freq, fftshift(abs(fft(pulse))), ...
    linspace(-freqlimit,freqlimit, N), ...
    fftshift(abs(fft(s_ambient(1:N)))), ':',...
    linspace(-freqlimit,freqlimit, length(time_corr)), ...
    fftshift(abs(fft(time_corr))), ':');
xlabel('[Hz]'); legend('|R|', 'ambient', 'FFT(corr)');
grid; axis([freq(1) freq(end) 1 1E5]);

%% Correlate using FFT and IFFT (fast correlate)
pulse_pad = [pulse zeros(size(pulse))];
FFT_pulse = fft(pulse_pad);
FFT_pulse_conj = conj(FFT_pulse);

FFT_pulse_conj_vDSPz = FFT_pulse_conj(1:M); % Just take half
FFT_pulse_conj_vDSPz(1) = FFT_pulse_conj(1) + 1j*FFT_pulse_conj(M+1);
% Write the flattened form (row vector) of FFT_pulse_conj_vDSPz to the header
FFT_pulse_conj_vDSPz_flattened = [ ...
    real(FFT_pulse_conj_vDSPz) imag(FFT_pulse_conj_vDSPz) ...
    ];
csvwrite(strcat(osxSrcDir, 'FFT_pulse_conj_vDSPz_flattened.h'), ...
    FFT_pulse_conj_vDSPz_flattened);

%unityDocDir = '~/github/OpenCVUnity/doc/';
%fname = strcat(unityDocDir, 'FFT_pulse_conj.mat');
%save(fname);

block_rms(1:N-M) = 0.0;

fft_corr = zeros(size(s_r) - [0 M]);
rx_pad = [zeros(size(pulse)) s_r(1:M)];
block_rms(1:M) = max(1, 0.5*sqrt(dot(rx_pad((M+1):end), rx_pad((M+1):end)))) ...
    * ones(1,M);
corr = ifft(fft(rx_pad) .* FFT_pulse_conj);
fft_corr(1:M) = corr(M + (1:M)); % Throw away the L half, save the R half
figure(3);
subplot(414); plot(t(1:length(fft_corr)), fft_corr);

for i=M:M:(N-2*M) % Feed subsequent received samples
    rx_pad((M+1):end) = s_r(i+(1:M));
    block_rms(i+(1:M)) = max(1, 0.5*sqrt(dot(rx_pad((M+1):end), ...
        rx_pad((M+1):end)))) ...
        * ones(1,M);
    corr = ifft(fft(rx_pad) .* FFT_pulse_conj);
    fft_corr(i + ((-M+1):M)) = fft_corr(i + ((-M+1):M)) + corr;
    subplot(414); plot(t(1:length(fft_corr)), fft_corr);
end
% Pick up the last block's L half (and add to the previous block's R half)
if isempty(i), i = M; else i = i + M; end
rx_pad((M+1):end) = s_r(i+(1:M));
corr = ifft(fft(rx_pad) .* FFT_pulse_conj);
fft_corr(i+((-M+1):0)) = fft_corr(i+ ((-M+1):0)) + corr(1:M);
err = time_corr - fft_corr;
if max(abs(err)) > 1E-9
    error('Fast correlation is wrong; e = %f', max(abs(err)))
    ha(4) = subplot(414); plot(t(1:length(err)), err);
    linkaxes(ha, 'x');
end

%block_rms = 0.5 * block_rms; % Set the threshold at half of RMS

ha(4) = subplot(414); plot(t(1:length(fft_corr)), fft_corr);
title('fast corr (overlapped)');
%linkaxes(ha, 'x');
%close(3);

%expensive_corr = ifft(fft(s_r) .* conj(fft(pulse)));
% if abs(txrx_corr - fft_corr) > 1E-9
%     error('Overlapped correlation is wrong!')
% end

%% Show the correlation
distance = t(1:length(fft_corr)) * c;
figure(2)
ha4(1) = subplot(311);
plot(distance, fft_corr, distance, block_rms);
title('FFT based correlation'); grid
xlabel ('Target relative distance [m]'); ylabel ('Correlation');

d2Corr = [0 diff(diff(fft_corr)) 0];
thresheld_d2Corr = max(-d2Corr, block_rms);

subplot(312);
plot(distance(1:length(d2Corr)), -d2Corr, ...
    ...distance(1:length(d2Corr)), -block_rms, ...
    distance(1:length(d2Corr)), thresheld_d2Corr);
grid;

med_corr = medfilt1(fft_corr, 3);
dCorr = fft_corr - med_corr;
thresheld_dCorr = max(dCorr, 0.5*block_rms);
subplot(313);
plot(distance(1:length(dCorr)), dCorr, ...
    ...distance(1:length(dCorr)), 0.5*block_rms, ...
    distance(1:length(dCorr)), thresheld_dCorr);
grid;

%close(2);
%close(3);

[~, k_dCor] = max(thresheld_dCorr);
[~, k_d2cor] = max(thresheld_d2Corr);
fprintf('diff of estimate = %d\n', k_dCor - k_d2cor);
d_est = Ts * k_dCor * c;
if abs(k_dCor - k_d2cor) <= 2
    fprintf('Have agreement on distance = %f\n', d_est);
else
    fprintf('Distance based on relative corr. d = %f\n', d_est);
end
return;

%% Analyze real data captured with Built-in mic
%fname = '~/radar/Qosx/from_self_macbook.csv';
fname = strcat(osxDataDir, 'from_ipod.csv');
%fname = '~/radar/QiOS/from_self_ipod.csv';
debug = csvread(fname); %load(fname);
heard = [debug(:,1); debug(:,2)];
figure(1);
%t = (1:length(heard)) / Fs;
subplot(211); plot(heard(:,1)); xlabel('[k]');
%ha(2) = subplot(212); plot(t, heard(:,2));
%linkaxes(ha, 'x');

% What is the frequency content of the heard signal?
freq = linspace(-freqlimit,freqlimit, length(heard));
FFT_heard = fft(heard);
subplot(212); plot(freq, abs(fftshift(FFT_heard)));
xlabel('[Hz]');

N = M * floor(length(heard) / M);

[s_ambient, Fs_a] = audioread('singing.wav', [1 N]);
s_r = heard(1:N,1)'; %+ s_ambient(1:N)'; % received = signal + noise + ambient

block_rms(1:N-M) = 0.0;

fft_corr = zeros(size(s_r) - [0 M]);
rx_pad = [zeros(size(pulse)) s_r(1:M)];
block_rms(1:M) = max(1, 0.5*sqrt(dot(rx_pad((M+1):end), rx_pad((M+1):end)))) ...
    * ones(1,M);
corr = ifft(fft(rx_pad) .* FFT_pulse_conj);
fft_corr(1:M) = corr(M + (1:M)); % Throw away the L half, save the R half
figure(3);
subplot(414); plot(t(1:length(fft_corr)), fft_corr);

for i=M:M:(N-2*M) % Feed subsequent received samples
    rx_pad((M+1):end) = s_r(i+(1:M));
    block_rms(i+(1:M)) = max(1, 0.5*sqrt(dot(rx_pad((M+1):end), ...
        rx_pad((M+1):end)))) ...
        * ones(1,M);
    corr = ifft(fft(rx_pad) .* FFT_pulse_conj);
    fft_corr(i + ((-M+1):M)) = fft_corr(i + ((-M+1):M)) + corr;
    subplot(414); plot(t(1:length(fft_corr)), fft_corr);
end
% Pick up the last block's L half (and add to the previous block's R half)
if isempty(i), i = M; else i = i + M; end
rx_pad((M+1):end) = s_r(i+(1:M));
corr = ifft(fft(rx_pad) .* FFT_pulse_conj);
fft_corr(i+((-M+1):0)) = fft_corr(i+ ((-M+1):0)) + corr(1:M);
close(3);

distance = t(1:length(fft_corr)) * c;
figure(4); clf
ha4(1) = subplot(311);
plot(distance, fft_corr, distance, block_rms);
title('FFT based correlation and rx signal power'); grid
xlabel ('Target relative distance [m]'); ylabel ('Correlation');

d2Corr = [0 diff(diff(fft_corr)) 0];
thresheld_d2Corr = max(-d2Corr, block_rms);

subplot(312);
plot(distance(1:length(d2Corr)), -d2Corr, ...
    ... distance(1:length(d2Corr)), -block_rms, ...
    distance(1:length(d2Corr)), thresheld_d2Corr);
title('THRESHOLD(d^2(corr)/dt)'); grid;

med_corr = medfilt1(fft_corr, 5);
dCorr = fft_corr - med_corr;
thresheld_dCorr = max(dCorr, 0.5*block_rms);

subplot(313);
plot(distance(1:length(dCorr)), dCorr, ...
    ... distance(1:length(d2Corr)), 0.5*block_rms, ...
    distance(1:length(d2Corr)), thresheld_dCorr);
title('THRESHOLD(corr - med5(corr))'); grid;

[~, k_dCor] = max(thresheld_d2Corr);
[~, k_d2cor] = max(thresheld_dCorr);
fprintf('diff of estimate = %d\n', k_dCor - k_d2cor);
d_est = Ts * k_dCor * c;
if abs(k_dCor - k_d2cor) <= 2
    fprintf('Have agreement on distance = %f\n', d_est);
else
    fprintf('Distance based on relative corr. d = %f\n', d_est);
end

return
