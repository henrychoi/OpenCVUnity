function matched_lfm
osxDataDir = '~/radar/Qosx/';
iosDataDir = '~/radar/QiOS/';
fname = strcat(iosDataDir, 'from_self_ipod.csv');
%fname = strcat(iosDataDir, 'from_macbook.csv');

AMBIENT_SOUND = 0
CHIRP_STRENGTH = 1

Fs = 44100;
Ts = 1/Fs;
% MacBook Pro fades in on the low end around 130 Hz
f0 = 5000; % lowest frequency of the chirp
b = 10000; % bandwidth
M = 2^10; % # samples to describe the original chirp (controls tau)
N = M * 2; %Total samples I need to collect
c = 335; % speed of sound [m/s]
taup = M * Ts;
Rmax = N * Ts * c;
freqlimit = 0.5 * Fs;

%fprintf('Detection range: %f m\n', 0.5 * Nsample * Ts * c);
fprintf('Range resolution: %f mm\n', 1000*c/b);

mu = b/taup; % compression ratio; @ t=taup, instantaneous F = f0 + b
%eps = 1.0e-16;

t = Ts * (0:N-1);
t_xmit = t(1:M);

% The ideal pulse, to estimate where the received waveform began
pulse = sin(2*pi * (f0 .* t_xmit + 0.5*mu .* t_xmit.^2) ...
    + 0.1*pi); % Start with non-zero amplitude
debug = csvread(fname); %load(fname);
heard = [debug(:,1); debug(:,2)];
% Correlate the heard against ideal
txrx_xcorr = xcorr(heard, pulse); % len = 2M+M+M - 1
[~, kcorr] = max(txrx_xcorr);
kcorr = kcorr - 2*M;

%kcorr = 745;
figure(1);
subplot(211); plot(1:N, heard); xlabel('[k]');
subplot(212);
if kcorr >= 0
    plot(1:M, heard(kcorr + (1:M)), 1:M, pulse * rms(heard(kcorr + (1:M))));
else
    plot(1:(M+kcorr), heard(1:(M+kcorr)), ...
        1:(M+kcorr), pulse((1-kcorr):M) * rms(heard(1:(M+kcorr))) );
end
xlabel('[k]');

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

