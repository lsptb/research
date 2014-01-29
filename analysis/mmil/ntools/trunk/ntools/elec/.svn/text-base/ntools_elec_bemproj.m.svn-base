function electrodes = ntools_elec_bemproj(ini_loc,ini_pos,dim,sph,bem_name,bem_surf_path)
% calculate the points on the bem surface

[temp_loc,temp_pos] = ntools_elec_223(ini_loc,ini_pos,dim(1),dim(2),sph);
elec = ntools_elec_grid(temp_loc,temp_pos,dim(1),dim(2));
l = bem_name(findstr('.',bem_name)+1:end);
if strcmp(l,'tri')
    [surf] = fs_read_trisurf([bem_surf_path bem_name]);
else 
    [surf] = fs_read_surf([bem_surf_path bem_name]);
end
ny.vert = surf.coords;
ny.tri = surf.faces;

brain = ny.vert;
if strcmp(sph,'lh')
    m = find(brain(:,1)>-5); % for left hemi
elseif strcmp(sph,'rh') 
    m = find(brain(:,1)<5); % for right hemi
end
brain(m,:) = [];
% x = elec(:,2);
y = elec(:,3);
z = elec(:,4);

for i = 1:length(y)
%     dst_x = abs(brain(:,1)-x(i));
    dst_y = abs(brain(:,2)-y(i));
    dst_z = abs(brain(:,3)-z(i));
    c = sqrt(dst_z.^2+dst_y.^2);
    index(i) = find(c==min(c));
end
electrodes = brain(index,:);

