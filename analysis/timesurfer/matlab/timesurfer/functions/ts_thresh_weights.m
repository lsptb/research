function [tw,tv]=thresh_weights(w,v,thresh,threshabsflag)
% thresh_weights.m: threshold weights
% [tw,tv] = thresh_weights(w,v,thresh,threshabsflag)
tw = w;
tv = v;
if(isempty(w))
  return;
end
try, threshabsflag; catch, threshabsflag = 1; end
if threshabsflag
  iw = find(abs(w)>thresh);
else
  iw = find(w>thresh);
end
tw = w(iw);
tv = v(iw);
return;
end
