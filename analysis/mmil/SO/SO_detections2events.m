function events = SO_detections2events(detections,times,labels)

events  = [];
for k   = 1:length(detections)
  if nargin > 2
    events(k).label = labels{k};
  else
    events(k).label = sprintf('channel %g',k);
  end
  events(k).time  = [times([detections(k).pospeak detections(k).negpeak])];
  events(k).type  = [1*ones(1,length(detections(k).pospeak)) 2*ones(1,length(detections(k).negpeak))];
end