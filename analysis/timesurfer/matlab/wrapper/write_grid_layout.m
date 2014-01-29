function layfile = write_grid_layout(filepath,sens,nr,nc,rpad,cpad)
%	write_grid_layout(layout1,data.sensor_info,chans,8,8,.01,.01)

badchan_idx     = find([sens.badchan]);
ignore_channels = [{sens(badchan_idx).label}];
max_group_size = 130;

fig = ts_iEEG_makefigs (sens,nr,nc,max_group_size,ignore_channels,0);

nr = []; nc = [];
for i = 1:length(fig)
	nr = max([fig(i).plots.row]);
	nc = max([fig(i).plots.column]);
	
	h = (.90 - nr*rpad) / nr;
	w = (.85 - nc*cpad) / nc;

	k = 1:nr;  y = [.90 - (k-1)*rpad - k*h]; 			% y pos
	k = 1:nc;  x = [.05 + (k-1)*cpad + (k-1)*w]; 	% x pos
	
	layfile{i} = sprintf('%s/layout_%s.lay',filepath,fig(i).name);
	fid = fopen(layfile{i},'wt');
	for j = 1:length(fig(i).plots)
		if any(badchan_idx==fig(i).plots(j).index), continue; end
		r 	= fig(i).plots(j).row;
		c 	= fig(i).plots(j).column;
		id 	= fig(i).plots(j).index;
		lbl = fig(i).plots(j).name;
	  fprintf(fid,'%g %.6f %.6f %.6f %.6f %s\n',id,x(c),y(r),w,h,lbl);% id x y width height label
	end
end
fclose(fid);
