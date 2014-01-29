setparallel(0,0);
pause(1);
ffp_start_record('/mnt/raid/adam','datalatency');
pause(1);
STATE = 0;
for i = 1:1000
  if STATE == 0
    Nspike_pretime(i) = gettime_nspike;
    Comedi_pretime(i) = gettime_comedi;
    setparallel(0,1);
    Nspike_posttime(i) = gettime_nspike;
    Comedi_posttime(i) = gettime_comedi;
    STATE = 1;
  elseif STATE == 1
    Nspike_pretime(i) = gettime_nspike;
    Comedi_pretime(i) = gettime_comedi;
    setparallel(0,0);
    Nspike_posttime(i) = gettime_nspike;
    Comedi_posttime(i) = gettime_comedi;
    STATE = 0;
  end
  pause(0.1);
end
ffp_stop_record;
