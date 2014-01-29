function [data,p] = gen_timecode(i,sr)

if nargin < 2 sr = 1e4; end

persistent waves;
persistent waves_sr;

    mult = 64;
    t       = linspace(0,1,sr+1);
    t = t(1:end-1);
    i = 2.*pi.*mult.*[32:-1:1];
    com = i'*t;
    waves = sin(com);
    waves_sr = sr;

    disp('entering gen_timecode')
    if nargin ==0
        i=0;
        while 1
            if i == 2^32-1; i = 0; end
            i=i+1;
            bin = (dec2bin(i,32) == '1');
            data = 0.99*bin*waves./sum(bin);
            p = audioplayer(data,sr);
            disp(['Sending timecode ' num2str(i)]);
            playblocking(p);
            pause(10);
            stop(p);
        end
    else
        data = modulate_timecode(i,sr);
        p = audioplayer(data,sr);
        play(p);
    end
