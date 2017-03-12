function matched_lfm
%% Simulate
tempDir = '~/Music/';
unityDocDir = '~/github/OpenCVUnity/doc/'
unityAssetsDir = '~/github/OpenCVUnity/Assets/'
soundDir = strcat(unityAssetsDir, 'CubeShipsFree/Sound/');
streamingAssetDir = strcat(unityAssetsDir, 'StreamingAssets/');
osxSrcDir = '~/github/OpenCVUnity/doc/CAInput/CH08_AUGraphInput/';

TRUE_DISTANCE = 7
AMBIENT_SOUND = 1
CHIRP_STRENGTH = 1
COMPLEX_SIGNAL = false
ROUND_TRIP = false
HIGH_PASS = false

if ROUND_TRIP, attenuation = 4; else attenuation = 2; end;
MULTIPATH_DISTANCE = TRUE_DISTANCE + 0.02;
Fs = 44100;
Ts = 1/Fs;
f0 = 50; % lowest frequency of the chirp
b = 22000; % 23 kHz BW (use all available BW)
M = 2^10; % # samples to describe the original chirp (controls tau)
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
winid = 4; % my own window
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
sigma_r = 0.1 * RMS_ambient;
A_u = CHIRP_STRENGTH * RMS_ambient * sqrt(2);

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

switch winid % determine proper window
    case 1
        win = hamming(M)';
    case 2
        win = kaiser(M, pi)';
    case 3
        win = chebwin(M, 60)';
    case 4
        win = ones(1,M);
        forte_d = round(200 * Fs / 48000);
        piano_d = round(50 * Fs / 48000);
        win(1:forte_d) = linspace(0,1, forte_d);
        win((M-piano_d+1):M) = linspace(1,0, piano_d);
end
%win = sqrt(win); % don't attenuate so much

% Form the transmitted signal
if COMPLEX_SIGNAL
    pulse = win .* exp(1j*2*pi * (f0 .* t_xmit + 0.5*mu .* t_xmit.^2));
else %only the real
    pulse = win .* cos(2*pi * (f0 .* t_xmit + 0.5*mu .* t_xmit.^2));
end

%% Write pulse to play back Unity
silence = zeros(size(pulse));
audiowrite(strcat(tempDir, 'lchirp.wav'), [pulse; silence]', Fs);
audiowrite(strcat(tempDir, 'rchirp.wav'), [silence; pulse]', Fs);

xmit = A_u * pulse; % Adjust pulse loudness
sound(xmit, Fs);% What does the chirp sound like?
%return;

k_true = round(Fs * scat_range(1) / c);
sdots = dot(xmit, pulse) / (scat_range(1).^2);
sdotn = zeros(1, 20);
for i=1:length(sdotn)
    sdotn(i) = dot(s_r((k_true + i - length(sdotn)/2) + (1:M)), pulse);
end
figure(4); subplot(212); plot(sdotn)

% Simulate multi-path
for j = 1:nscat
    range = scat_range(j);
    % Discretize the distance into the sample instance
    if ROUND_TRIP, range = 2*range; end;
    
    t_R = range/c; % nominal wave travel time
    x(j, (1:M) + round(t_R * Fs)) = xmit;
    
    % scale receoved signal and accumulate
    s_r = s_r + (scat_rcs(j) / (range)^attenuation) * x(j,:);
end
sdotn = dot(s_r(k_true + (1:M)), pulse);

figure(2);
subplot(211); plot(t_xmit, xmit, t(1:N), s_ambient((1:N)), ':', t(1:N), rcv_noise, '--');
title('truth'); legend('xmit', 'ambient', 'noise');
subplot(212); semilogy(t(1:M), abs(xmit) .^2, t,abs(s_r).^2, ':'); %grid
title('Sound amplitude'); legend('xmit', 'received'); xlabel('[s]');
axis([-inf inf 1E-3, inf])

% What do correlations look like?
%return

if HIGH_PASS %% What does the high pass of received signal look like?
    fname = '~/Documents/MATLAB/Hhp.mat';
    %save(fname, 'Hhp');
    load(fname);
    figure(1);
    subplot(211); plot(Hhp.Numerator);
    % Fast convolution, using padded filter kernel
    hpf = zeros(1,2*M);%zero pad
    %hpf(1:(Hhp.order/2+1)) = Hhp.Numerator((Hhp.order/2+1):end);
    %hpf((end-Hhp.order/2+1):end) = Hhp.Numerator(1:(Hhp.order/2));
    %hpf(M + (-Hhp.order/2:Hhp.order/2)) = Hhp.Numerator;
    hpf(1:length(Hhp.Numerator)) = Hhp.Numerator;
    %hpf = circshift(hpf, -45);
    FFT_filter = fft(hpf);
    subplot(211);
    freq_pad = linspace(-freqlimit,freqlimit,2*M);
    plot(freq_pad, abs(fftshift(FFT_filter)));
    % fast convolution
    s_filt = zeros(size(s_r));
    rx_pad = [s_r(1:M) zeros(1,M)];
    subplot(211); plot(abs(fft(rx_pad)));
    Fresp = fft(rx_pad) .* FFT_filter; corr = ifft(Fresp);
    s_filt(1:2*M) = corr;

    figure(3);
    subplot(411); plot(filtfilt(Hhp.Numerator, 1, s_r));
    subplot(412); plot(freq_pad, fftshift(Fresp));
    subplot(413); plot(t(1:length(corr)), corr);
    subplot(414); plot(t(1:length(s_filt)), s_filt);

    for i=M:M:(N-2*M) % Feed subsequent received samples
        rx_pad(1:M) = s_r(i+(1:M));
        Fresp = fft(rx_pad) .* FFT_filter; corr = ifft(Fresp);
        s_filt(i + (1:2*M)) = s_filt(i + (1:2*M)) + corr;
        subplot(412); plot(freq_pad, Fresp);
        subplot(413); plot(t(1:length(corr)), corr);
        subplot(414); plot(t(1:length(s_filt)), s_filt);
    end
    % Pick up the last block's L half (and add to the previous block's R half)
    i = i + M;
    rx_pad(1:M) = s_r(i+(1:M));
    Fresp = fft(rx_pad) .* FFT_filter; corr = ifft(Fresp);
    s_filt(i + (1:M)) = s_filt(i + (1:M)) + corr(1:M);
    subplot(412); plot(freq_pad, Fresp);
    subplot(413); plot(t(1:length(corr)), corr);
    subplot(414); plot(t(1:length(s_filt)), s_filt);
end %HIGH_PASS

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
save(strcat(unityDocDir, 'FFT_pulse_conj.mat'), 'FFT_pulse_conj');
% Using jsonlab-1.5, downloaded from Matlab file exchange
savejson('', FFT_pulse_conj_vDSPz, ...
    strcat(streamingAssetDir, 'FFT_pulse_conj_vDSPz.json'));
% Write the flattened form (row vector) of FFT_pulse_conj_vDSPz to the header
FFT_pulse_conj_vDSPz_flattened = [ ...
    real(FFT_pulse_conj_vDSPz) imag(FFT_pulse_conj_vDSPz) ...
    ];
csvwrite(strcat(osxSrcDir, 'FFT_pulse_conj_vDSPz_flattened.h'), ...
    FFT_pulse_conj_vDSPz_flattened);

fft_corr = zeros(size(s_r) - [0 M]);
rx_pad = [zeros(size(pulse)) s_r(1:M)];
corr = ifft(fft(rx_pad) .* FFT_pulse_conj);
fft_corr(1:M) = corr(M + (1:M)); % Throw away the L half, save the R half
figure(3);
subplot(414); plot(t(1:length(fft_corr)), fft_corr);

for i=M:M:(N-2*M) % Feed subsequent received samples
    rx_pad((M+1):end) = s_r(i+(1:M));
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
if ROUND_TRIP, distance = distance/2; end;

figure(4)
ha4(1) = subplot(311);
plot(distance, fft_corr); title('FFT based correlation'); grid
xlabel ('Target relative distance [m]'); ylabel ('Correlation');

d2Corr = [0 diff(diff(fft_corr)) 0];
subplot(312); plot(distance(1:length(d2Corr)), d2Corr);
grid;

% What does the smoothed signal look?

if HIGH_PASS % The filtered version
    s_filt = zeros(size(d2Corr));
    rx_pad = [d2Corr(1:M) zeros(1,M)];

    s_filt = zeros(size(fft_corr));
    rx_pad = [fft_corr(1:M) zeros(1,M)];
    Fresp = fft(rx_pad) .* FFT_filter; corr = ifft(Fresp);
    s_filt(1:2*M) = corr;
    figure(3);
    subplot(411); plot(filtfilt(Hhp.Numerator, 1, fft_corr));
    subplot(412); plot(freq_pad, fftshift(Fresp));
    subplot(413); plot(t(1:length(corr)), corr);
    subplot(414); plot(t(1:length(s_filt)), s_filt);

    for i=M:M:(N-3*M) % Feed subsequent received samples
        rx_pad(1:M) = d2Corr(i+(1:M));
        Fresp = fft(rx_pad) .* FFT_filter; corr = ifft(Fresp);
        s_filt(i + (1:2*M)) = s_filt(i + (1:2*M)) + corr;
        subplot(412); plot(freq_pad, Fresp);
        subplot(413); plot(t(1:length(corr)), corr);
        subplot(414); plot(t(1:length(s_filt)), s_filt);
    end
    % Pick up the last block's L half (and add to the previous block's R half)
    if isempty(i), i = M; else i = i + M; end
    rx_pad(1:M) = d2Corr(i+(1:M));
    Fresp = fft(rx_pad) .* FFT_filter; corr = ifft(Fresp);
    s_filt(i + (1:M)) = s_filt(i + (1:M)) + corr(1:M);
    subplot(412); plot(freq_pad, Fresp);
    subplot(413); plot(t(1:length(corr)), corr);
    subplot(414); plot(t(1:length(s_filt)), s_filt);

    figure(4); subplot(313); plot(distance, s_filt); title('FFT correlation filtered');
    grid
end
close(2);
close(3);

med_corr = medfilt1(fft_corr, 5);
dCorr = fft_corr - med_corr;
subplot(313); plot(distance(1:length(dCorr)), dCorr);
grid;


[~, k_dCor] = max(dCorr);
[~, k_d2cor] = min(d2Corr);
fprintf('diff of estimate = %d\n', k_dCor - k_d2cor);
d_est = Ts * k_dCor * c;
if abs(k_dCor - k_d2cor) <= 2
    fprintf('Have agreement on distance = %f\n', d_est);
else
    fprintf('Distance based on relative corr. d = %f\n', d_est);
end

%% Analyze real data captured with Built-in mic
fname = '~/github/OpenCVUnity/doc/heard.mat';
load(fname);
figure(1);
t = (1:length(heard)) / Fs;
ha(1) = subplot(211); plot(t, heard(:,1)); xlabel('[s]');
%ha(2) = subplot(212); plot(t, heard(:,2));
%linkaxes(ha, 'x');

% What is the frequency content of the heard signal?
freq = linspace(-freqlimit,freqlimit, length(heard));
FFT_heard = fft(heard);
subplot(212); plot(freq, abs(fftshift(FFT_heard)));
xlabel('[Hz]');

N = M * floor(length(heard) / M);
[s_ambient, Fs_a] = audioread('singing.wav', [1 N]);
s_r = heard(1:N,1)' + s_ambient(1:N)';
fft_corr = zeros(size(s_r) - [0 M]);
rx_pad = [zeros(size(pulse)) s_r(1:M)];
corr = ifft(fft(rx_pad) .* FFT_pulse_conj);
fft_corr(1:M) = corr(M + (1:M)); % Throw away the L half, save the R half
figure(3);
subplot(414); plot(t(1:length(fft_corr)), fft_corr);

for i=M:M:(N-2*M) % Feed subsequent received samples
    rx_pad((M+1):end) = s_r(i+(1:M));
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
plot(distance, fft_corr); title('FFT based correlation'); grid
xlabel ('Target relative distance [m]'); ylabel ('Correlation');

d2Corr = [0 diff(diff(fft_corr)) 0];
subplot(312); plot(distance(1:length(d2Corr)), d2Corr);
grid;

med_corr = medfilt1(fft_corr, 5);
dCorr = fft_corr - med_corr;
subplot(313); plot(distance(1:length(dCorr)), dCorr);
grid;

[~, k_dCor] = max(dCorr);
[~, k_d2cor] = min(d2Corr);
fprintf('diff of estimate = %d\n', k_dCor - k_d2cor);
d_est = Ts * k_dCor * c;
if abs(k_dCor - k_d2cor) <= 2
    fprintf('Have agreement on distance = %f\n', d_est);
else
    fprintf('Distance based on relative corr. d = %f\n', d_est);
end

return
