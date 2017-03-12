%% Capture audio of a chirp
for i=0:(audiodevinfo(1)-1)
    if isempty(strfind(audiodevinfo(1,i), 'Built'))
        continue;
    end
    break;
end
Fs = 44100
% OK, now I have the ID of the recorder
recorder = audiorecorder(Fs, 16, 2, i)
disp('recording...');
recordblocking(recorder, 5);
disp('recording done');

play(recorder);
heard = getaudiodata(recorder);
plot(heard);

fname = '~/github/OpenCVUnity/doc/heard.mat';
save(fname, 'heard');
