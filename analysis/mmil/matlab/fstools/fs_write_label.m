function ok = write_label(lindex, lxyz, lvals, labelfile, subjid)
% ok = write_label(lindex, lxzy, lvals, labelfile, <subjid>)
%
% vertex numbers in lindex should be 1-based (matlab style)
%
ok = 0;

if(nargin ~= 4 & nargin ~= 5)
  fprintf('ok = write_label(lindex, lxzy, lvals, labelfile, <subjid>)\n');
  return;
end

if(exist('subid') ~= 1) subjid = ''; end

if(isempty(lindex) & isempty(lxyz))
  fprintf('ERROR: both lindex and lxyz are empty.\n');
  return;
end

if(~isempty(lindex) & ~isempty(lxyz))
  npoints1 = length(lindex); 
  npoints2 = size(lxyz,1); 
  if(npoints1 ~= npoints2)
    fprintf('ERROR: lindex and lxyz have different lengths.\n');
    return;
  end
  npoints = npoints1;
elseif(~isempty(lindex))
  npoints = length(lindex); 
  lxyz = zeros(npoints,3); 
elseif(~isempty(lxyz))
  npoints = length(lxyz); 
  lindex = zeros(npoints,1); 
end

if(size(lxyz,2) ~= 3)
  fprintf('ERROR: lxyz does not have 3 columns\n');
  return;
end

if(~isempty(lvals)) 
  if(npoints ~= length(lvals))
    fprintf('ERROR: length of lvals inconsistent\n');
    return;
  end
else
  lvals  = zeros(npoints,1); 
end

% open as an ascii file
fid = fopen(labelfile, 'w') ;
if(fid == -1)
  fprintf('ERROR: could not open %s\n',labelfile);
  return;
end

fprintf(fid,'#!ascii label, from subject %s \n',subjid);
fprintf(fid,'%d\n',npoints);

% Make sure they are npoints by 1 %
lindex = reshape(lindex,[npoints 1]) - 1; % change to zero base
lxyz   = reshape(lxyz,[npoints 3]);
lvals  = reshape(lvals,[npoints 1]);

l = [lindex lxyz lvals];
fprintf(fid,'%d %f %f %f %f\n',l') ;

fclose(fid) ;

ok = 1;

return;
