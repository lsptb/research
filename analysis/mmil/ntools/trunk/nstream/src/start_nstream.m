nstream_add_auxdsp('10.1.2.11');
nstream_add_auxdsp('10.1.2.12');
nstream_add_auxdsp('10.1.2.13');
nstream_add_auxdsp('10.1.2.14');
nstream_add_auxdsp('10.1.2.15');
nstream_add_auxdsp('10.1.2.16');
nstream_add_auxdsp('10.1.2.17');
nstream_add_auxdsp('10.1.2.18');
nstream_add_auxdsp('10.1.2.21');
nstream_add_auxdsp('10.1.2.22');
nstream_add_auxdsp('10.1.2.23');
nstream_add_auxdsp('10.1.2.24');
nstream_add_auxdsp('10.1.2.25');
nstream_add_auxdsp('10.1.2.26');
nstream_add_auxdsp('10.1.2.27');
nstream_add_auxdsp('10.1.2.28');
i=1;

for twice = 0:1
    for group = 0:3
        for division = 0:7
            for channel = 0:3
                hw_num = (division * 16) + (group * 4) + channel;
                nstream_set_channel(i, hw_num);
                i=i+1;
            end
        end
    end
end            

nstream_start_acquire

