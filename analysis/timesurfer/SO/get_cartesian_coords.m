% purpose: get X,Y coords in 2D plane for sensor data
% written by JSS on 07-Oct-2010

function [X,Y] = get_cartesian_coords(data,layout)
if nargin < 2, layout = []; end
nchan   = data.num_sensors;
pos     = zeros(nchan,3); % (x,y,z) for each channel
if isempty(layout)
  T     = data.coor_trans.device2head;
  for k = 1:nchan
    loc         = data.sensor_info(k).loc;
    if any(strmatch('grad',data.sensor_info(k).typestring)) || ...
       any(strmatch('mag' ,data.sensor_info(k).typestring))
      loc       = T*loc;
    end  
    pos(k,1:3)  = loc(1:3,4);
  end
  method = 'stereographic';%params.method; % gnomic, stereographic, ortographic, inverse, polar
  prj    = elproj(pos, method); % * [0 1; -1 0];
            % ELPROJ makes a azimuthal projection of a 3D electrode cloud
            %  on a plane tangent to the sphere fitted through the electrodes
            %  the projection is along the z-axis
  X = prj(:,1);   % x-coordinates
  Y = prj(:,2);   % y-coordinates  
elseif ischar(layout) && exist(layout,'file')
  [chNum,X,Y,Width,Height,Lbl,Rem] = textread(layout,'%f %f %f %f %f %q %q');
  for i=1:length(Lbl)
    if ~isempty(Rem{i})
      % this ensures that channel names with a space in them are also supported (i.e. Neuromag)
      Lbl{i} = [Lbl{i} ' ' Rem{i}];
    end
  end
  [sel,jnk] = match_str(Lbl,{data.sensor_info.label});
  X = X(sel);
  Y = Y(sel);
%   lay.pos    = [X Y];
%   lay.width  = Width;
%   lay.height = Height;
%   lay.label  = Lbl;
elseif isstruct(layout)
  X = lay.pos(:,1);
  Y = lay.pos(:,2);
else
  error('Layout not found.');
end