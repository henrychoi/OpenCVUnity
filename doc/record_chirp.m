%% Capture audio of a chirp
for i=0:(audiodevinfo(1)-1)
    if isempty(strfind(audiodevinfo(1,i), 'Built'))
        continue;
    end
    break;
end
% OK, now I have the ID of the recorder
recorder = audiorecorder(Fs, 16, 2, i)
disp('recording...');
recordblocking(recorder, 5);
disp('recording done');

play(recorder);
heard = getaudiodata(recorder);
plot(heard);
figure(1);
subplot(211); plot(heard(:,1));
subplot(212); plot(heard(:,2));

fname = '~/radar/txrx.mat';
save(fname, 'FFT_xmit_conj', 'xmit', 'heard');
