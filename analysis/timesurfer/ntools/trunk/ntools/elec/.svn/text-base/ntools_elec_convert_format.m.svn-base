function [nFile del_flag]= ntools_elec_convert_format(path,file)
% --------------------- convert the format
% if the del_flag is 1, then delete the converted file after running


file_noext = file(1:findstr(file,'.')-1);
file_ext = file(findstr(file,'.')+1:end);
if strcmp(file_ext,'nii.gz')
    gunzip([path file]);
    nFile = [file_noext '.nii'];
    del_flag = 1;
%     delete([path file]);
elseif strcmp(file_ext,'nii') || strcmp(file_ext,'img')
    nFile = file;
    del_flag = 0;
else
    convert = sprintf('mri_convert -ot nii %s%s %s%s.nii',path,file,path,file_noext);
    system(convert);
    nFile = [file_noext,'.nii'];
    del_flag = 1;
%     delete([path file]);
end
