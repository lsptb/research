%
% modulate_timecode(number, sr)
%
% encodes the supplied 32 bit integer in one second
% of analog data.
% 
% ex:
%
% soundsc(modulate_timecode(2^32-1, 44100), 44100, 16);
%
% generates one second of 16 bit audio data at the
% given sampling rate.  the first 900ms of the data 
% consists of a tone which is comprised of 1-32 
% sinewaves at 1-32 different frequencies representing
% 1-32 different bits.  the presence or absence of a
% given sinewave indicates the state of a given
% digital bit.
% 
function data=modulate_timecode(integer, sr)
% spacing in Hz between indicator frequencies
mult = 64;
%sr = 44100;
persistent waves;
persistent waves_sr;

if (isempty(waves_sr) || waves_sr ~= sr)
    disp('creating waves');
    t       = linspace(0,1,sr+1);
    t = t(1:end-1);
    i = 2.*pi.*mult.*[32:-1:1];
    com = i'*t;
    waves = sin(com);
    waves_sr = sr;
end

bin = (dec2bin(integer,32) == '1');
data = bin*waves;
data = data./(max(abs(data))+1);

% data(round(sr.*.9):sr) = 0;

end


