
h = plot(zeros(1,1e4));
axis([0 1e4 -65536 65536]);

while 1
%data = get_analog_state(10000);
disp('to:')
t = gettime_nspike
disp('from:')
t2 = t - 300

data = getrawint_nspike(t2, t);
%data = getraw_comedi_all;
set(h, 'YData', data(1,:));
axis([0 1e4 -65536 65536]);
pause(0.01);
end

