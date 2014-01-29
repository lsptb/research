function events = detections2events(detections)
t       = detections(1).tstart:1/detections(1).sfreq:detections(1).tstop;
events  = [];
for k   = 1:length(detections)
  events(k).label   = detections(k).label;
  events(k).time    = [t([detections(k).pospeak detections(k).negpeak])];
  events(k).latency = [detections(k).pospeak detections(k).negpeak];
  events(k).type    = [1*ones(1,length(detections(k).pospeak)) 2*ones(1,length(detections(k).negpeak))];
end
