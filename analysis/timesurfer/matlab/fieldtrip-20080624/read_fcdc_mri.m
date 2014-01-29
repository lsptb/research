function [mri] = read_fcdc_mri(filename);

% READ_FCDC_MRI reads anatomical and functional MRI data from different
% file formats. The output data is structured in such a way that it is
% comparable to a FieldTrip source reconstruction.
%
% Use as
%   [mri] = read_fcdc_mri(filename)
%
% The output MRI may have a homogenous transformation matrix that converts
% the coordinates of each voxel (in xgrid/ygrid/zgrid) into head
% coordinates. 
%
% See also READ_FCDC_DATA, READ_FCDC_HEADER, READ_FCDC_EVENT

% Copyright (C) 2004, Robert Oostenveld
%
% $Log: read_fcdc_mri.m,v $
% Revision 1.16  2007/05/06 09:09:18  roboos
% added support for nifti, requires SPM5
%
% Revision 1.15  2007/02/21 10:00:37  roboos
% give a more meaningfull warning when trying to read a mnc file without spm2 on the path
%
% Revision 1.14  2006/07/27 08:01:32  roboos
% use the hastoolbox function to detect the required toolboxes
%
% Revision 1.13  2006/04/19 15:46:43  roboos
% fixed bug for dicom when files are in another directory
%
% Revision 1.12  2006/04/19 07:54:51  roboos
% added support foir multifile DICOM format (searching for similar files)
% removed the trivial x/y/zgrid in the output structure
%
% Revision 1.11  2006/01/05 13:27:59  roboos
% change for avw_img_read(): force reading LAS* and flip the 1st dimension
%
% Revision 1.10  2005/09/08 07:32:48  roboos
% changed & into &&
% modified hasafni to check for BrikInfo instead of WriteBrik
%
% Revision 1.9  2005/08/17 19:35:30  roboos
% added a check for the presence of the file, give error if not present
%
% Revision 1.8  2005/07/20 15:33:35  roboos
% added a homogenous transformation matrix to AFNI import section
% I am not sure whether it is correct, but it works like this for the TTatlas
%
% Revision 1.7  2005/05/17 17:50:38  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.6  2004/10/28 09:41:20  roboos
% added preliminary support for AFNI mri files (no transformation matrix yet) using external afni_matlab toolbox
%
% Revision 1.5  2004/09/06 13:29:05  roboos
% fixed bug in transform for neuromag
%
% Revision 1.4  2004/09/03 09:17:49  roboos
% added support for neuromag MRI from fif file
%
% Revision 1.3  2004/08/26 12:12:50  roboos
% added transformation matrix for ASA files
%
% Revision 1.2  2004/08/26 11:57:22  roboos
% added support for MINC using SPM
% added support for Analyze using SPM
%
% Revision 1.1  2004/08/24 13:52:50  roboos
% new implementation, currently implemented are CTF, ASA and Analyze (with mri-toolbox)
%

% test for the presence of some external functions from other toolboxes
hasmri  = hastoolbox('mri');     % from Darren Weber, see http://eeg.sourceforge.net/
hasspm2 = hastoolbox('spm2');    % see http://www.fil.ion.ucl.ac.uk/spm/
hasspm5 = hastoolbox('spm5');    % see http://www.fil.ion.ucl.ac.uk/spm/
hasspm = (hasspm2 || hasspm5);
hasafni = hastoolbox('afni');    % see http://afni.nimh.nih.gov/

% test whether the file exists
if ~exist(filename)
  error(sprintf('file ''%s'' does not exist', filename));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if filetype(filename, 'ctf_mri')
  [img, hdr] = read_ctf_mri(filename);
  transform = hdr.transformMRI2Head;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif filetype(filename, 'asa_mri')
  [img, seg, hdr] = read_asa_mri(filename);
  transform = hdr.transformMRI2Head;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif filetype(filename, 'minc') 
  if ~hasspm
    error('the SPM2 or SPM5 toolbox is required to read *.mnc files');
  end
  % use the functions from SPM
  hdr = spm_vol_minc(filename);
  img = spm_read_vols(hdr);
  transform = hdr.mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif filetype(filename, 'nifti') 
  if ~hasspm5
    error('the SPM5 toolbox is required to read *.nii files');
  end
  % use the functions from SPM
  hdr = spm_vol_nifti(filename);
  img = spm_read_vols(hdr);
  transform = hdr.mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (filetype(filename, 'analyze_img') || filetype(filename, 'analyze_hdr')) && hasspm
  % use the image file instead of the header
  filename((end-2):end) = 'img';
  % use the functions from SPM to read the Analyze MRI
  hdr = spm_vol(filename);
  img = spm_read_vols(hdr);
  transform = hdr.mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (filetype(filename, 'analyze_hdr') || filetype(filename, 'analyze_img')) && hasmri
  % use the functions from Darren Weber's mri_toolbox to read the Analyze MRI
  avw = avw_img_read(filename, 0); % returned volume is LAS*
  img = avw.img;
  hdr = avw.hdr;
  % The default Analyze orientation is axial unflipped (LAS*), which means
  % that the resulting volume is according to the radiological convention.
  % Most other fMRI and EEG/MEG software (except Mayo/Analyze) uses
  % neurological conventions and a right-handed coordinate system, hence
  % the first axis of the 3D volume (right-left) should be flipped to make
  % the coordinate system comparable to SPM
  warning('flipping 1st dimension (L-R) to obtain volume in neurological convention');
  img = flipdim(img, 1);
  % FIXME: here I should also implement a homogenous transformation matrix,
  % using the voxel dimensions that are specified in hdr.dime.pixdim

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (filetype(filename, 'afni_brik') || filetype(filename, 'afni_head')) && hasafni
  [err, img, hdr, ErrMessage] = BrikLoad(filename);
  if err
    error('could not read AFNI file');
  end

  % FIXME: this should be checked, but I only have a single BRIK file
  % construct the homogenous transformation matrix that defines the axes
  warning('homogenous transformation might be incorrect for AFNI file');
  transform        = eye(4);
  transform(1:3,4) = hdr.ORIGIN(:);
  transform(1,1)   = hdr.DELTA(1);
  transform(2,2)   = hdr.DELTA(2);
  transform(3,3)   = hdr.DELTA(3);
  
  % FIXME: I am not sure about the "RAI" image orientation
  img = flipdim(img,1);
  img = flipdim(img,2);
  dim = size(img);
  transform(1,4) = -dim(1) - transform(1,4);
  transform(2,4) = -dim(2) - transform(2,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif filetype(filename, 'neuromag_fif')
  % use the neuromag functions to read the Neuromag MRI
  % FIXME this needs to be checked to ensure consistency with the FieldTrip definition of volume data
  [img,coords]  =loadmri(filename);
  dev = loadtrans(filename,'MRI','HEAD');
  transform = dev*coords;
  hdr.coords = coords;
  hdr.dev = dev;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif filetype(filename, 'dicom')
  % this uses the Image processing toolbox
  % the DICOM file probably represents a stack of slices
  [p, f] = fileparts(filename);
  tok = tokenize(f, '.');
  for i=5:length(tok)
    tok{i} = '*';  % this works for the Siemens scanners at the FCDC
  end
  filename = sprintf('%s.', tok{:});  % reconstruct the filename with wildcards and '.' between the segments
  filename = filename(1:end-1);       % remove the last '.'
  dirlist  = dir(fullfile(p, filename));
  dirlist  = {dirlist.name};
  for i=1:length(dirlist)
    filename = char(fullfile(p, dirlist{i}));
    fprintf('reading ''%s''\n', filename);
    info       = dicominfo(filename);
    img(:,:,i) = dicomread(info);
    hdr(i)     = info;
    if i==1
      % this pre-allocates enough space for the other slices
      img(1,1,length(dirlist)) = 0;
    end
  end
  % reorder the slices
  [z, indx]   = sort(cell2mat({hdr.SliceLocation}));
  hdr = hdr(indx);
  img = img(:,:,indx);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  error(sprintf('unrecognized filetype of ''%s''', filename));
end

% set up the axes of the volume in voxel coordinates
nx = size(img,1);
ny = size(img,2);
nz = size(img,3);
mri.dim = [nx ny nz];
% store the anatomical data
mri.anatomy = img;
% store the header with all fileformat specific details
mri.hdr = hdr;
try
  % if present, store the homogenous transformation matrix
  mri.transform = transform;
end

