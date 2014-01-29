if exist('sens','var')
elseif exist('data','var'),           sens = data.sensor_info;
elseif exist('epoch_data','var'),     sens = epoch_data.sensor_info;
elseif exist('avg_data','var'),       sens = avg_data.sensor_info;
elseif exist('timefreq_data','var'),  sens = timefreq_data.sensor_info;
elseif exist('plv_data','var'),       sens = plv_data.sensor_info;
elseif exist('flip','var') && isfield(flip,'sensor_info'), sens = flip(1).sensor_info;
end
if exist('sens','var')
  if ~exist('zshift','var'), zshift = .08; end
  if     exist('highlight','var') && isnumeric(highlight)
    highlight = highlight(highlight<=length(sens));
  elseif exist('highlight','var') && iscell(highlight)
    [highlight,jnk] = match_str({sens.label},highlight);
  elseif exist('highlight','var') && ischar(highlight)
    highlight = strmatch(highlight,{sens.label});
  else
    highlight = 1:length(sens);
  end
%   figure;
  x = arrayfun(@(x)(x.loc(1,4)),sens,'UniformOutput',true)';
  y = arrayfun(@(x)(x.loc(2,4)),sens,'UniformOutput',true)';
  z = arrayfun(@(x)(x.loc(3,4)),sens,'UniformOutput',true)';
  z = z - zshift;
  h = plot3(x,y,z,'b.'); set(h,'MarkerSize',4,'LineWidth',.4);
  for i=1:length(highlight)
    if ~exist('hlcolor','var'), hlcolor='k'; end
    htext = text(x(highlight(i)),y(highlight(i)),z(highlight(i)),'o','FontSize',8,'Color',hlcolor);
%     htext = text(x(highlight(i)),y(highlight(i)),z(highlight(i)),sens(highlight(i)).label,'FontSize',7,'Color',hlcolor);
  end
  axis([-.2 .2 -.2 .2 -.2 .1])
  grid on
  vline(0,'g')
  hline(0,'r')
end