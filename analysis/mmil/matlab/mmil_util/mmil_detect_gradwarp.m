function unwarpflag = mmil_detect_gradwarp(dcmfname)
%function unwarpflag = mmil_detect_gradwarp(dcmfname)
%
% loads dicom file, looks for zeros on edges (indicates gradwarp applied)
%   unwarpflag=1  (gradwarp has been applied)
%   unwarpflag=0  (gradwarp has not been applied)
%

thresh = 0.1;

if ~exist(dcmfname,'file')
  error('file %s not found');
end;
try
  im = dicomread(dcmfname);
catch
  error('unable to read file %s as dicom',dcmfname);
end;

border = [rowvec(im(1,:)),rowvec(im(end,:)),...
          rowvec(im(2:end-1,1)),rowvec(im(2:end-1,end))];

if length(find(border==0))/length(border) > thresh
  unwarpflag = 1;
else
  unwarpflag = 0;
end;

