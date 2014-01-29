run('exp_defn_adam');

last_time = nstream_gettime_nspike;

while 1

    dio_events = nstream_getdio_ts(last_time);
    pause(0.1);
    if (size(dio_events,2) > 0)
        last_time = dio_events(1,end);
        pause(2);
        for i = 1:size(dio_events,2)
            timestamp = dio_events(1,i);
            value = dio_events(4,i);
            disp(['dio event: timestamp=' num2str(dio_events(1,i)) ' value=' num2str(dio_events(4,i))]);
            eventdata = nstream_getrawint_nspike(timestamp - 1500, timestamp + 300);
            plot(eventdata(1,:));
        end
    end
end
