

STATE_LINE_LOW  = 0;
STATE_LINE_HIGH = 1;
state = STATE_LINE_LOW;

setparallel(0,0);
sfid = openserial;
pause(1);
ffp_start_record('/mnt/raid/adam','displaysensor_nspike_load')
pause(1);
%N = 3e4 - 1;
%h = plot(zeros(1,N));
%axis([0 N 0 65536])
x = [];
N = 1;

while 1
    d = getraw_nspike(N);
    s = d(1,1);
    % x = [x s];
%set(h,'YData',d(1,:));
    if s > 20000 & state == STATE_LINE_LOW
        disp('state change: high');
        setparallel(0,1);
        writeserial(sfid,'oId=2;oStatus=1;#');
        state = STATE_LINE_HIGH;
    elseif s < -20000 & state == STATE_LINE_HIGH
        disp('state change: low');
        setparallel(0,0);
        writeserial(sfid,'oId=2;oStatus=0;#');
        state = STATE_LINE_LOW;
    %elseif s < 6e4 & state == STATE_LINE_LOW
        %disp('No state change.  State is low')
    %elseif s > 6e4 & state == STATE_LINE_HIGH
        %disp('No state change.  State is high');
    end
    %pause(.01);
end

ffp_stop_record;
