function ini_name = ntools_elec_get_elec_name(ini_cell)
% suppose the elec was named as Letter+Number, this program just take the
% letter out
name = regexp(ini_cell(:,1),'[A-Za-z]*[^\d*]','match');
for i=1:length(name)
    ini_name(i) = name{i};
end
ini_name = unique(ini_name);