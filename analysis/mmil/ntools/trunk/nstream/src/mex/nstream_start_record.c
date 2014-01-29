#include <stdio.h>
#include <stdlib.h>
/*----------------------------------------------------------------------*/
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <string.h>
#include "nstreamclient.h"
#include "nstream.h"
#include "mex.h"

#define NAME_LENGTH 256
/*----------------------------------------------------------------------
% Y = nstream_start_record(recording_path, recording_filename_root, 
%	decimation_factor, comedi_num_channels_to_write, 
%	nspike_num_channels_to_write_low, nspike_num_channels_to_write_high,
%	nspike_high_channels_offset)
%
% nstream_start_record starts and stops recording 
%
% Inputs:	 recording_path          = directory in which to place this recording 
%		     recording_filename_root = root for generating filenames
%            decimation_factor       = (optional) decimate by a given integer on save
%	     comedi_num_channels     = (optional)  number of comedi channels to write
%
% Example:  nstream_start_record('/mnt/raid/Adam/010108/001', 'rec001')
%           creates: 
%                     /mnt/raid/Adam/010108/001/rec001.nspike.dat 
%                     /mnt/raid/Adam/010108/001/rec001.comedi.dat 
%                     /mnt/raid/Adam/010108/001/rec001.dio.txt 
%
%   Notes:  decimation factor defaults to 1, which is full rate 30k
%           if you enable decimation, you must configure filters on
%           the channels such that spectral content above the new
%           effective sampling rate after decimation is filtered out
%           otherwise there will be alaising in the saved data.
%
%	Defaults:
%	    comedi_num_channels_to_write 	8
%	    nspike_num_channels_to_write_low	256
%	    nspike_num_channels_to_write_high	0
%	    nspike_high_channels_offset 	0
%
% Outputs:  Y  =  
%
 ----------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char recording_path[NAME_LENGTH];
  char recording_filename_root[NAME_LENGTH];
  char error_text[NSTREAMCLIENT_ERROR_LENGTH];
  int decimation_factor = 1;
  int comedi_num_channels_to_write = NSTREAM_NUM_COMEDI_CHANNELS;
  int nspike_num_channels_to_write_low = NSTREAM_NUM_NSPIKE_CHANNELS;
  int nspike_num_channels_to_write_high = 0;
  int nspike_high_channels_offset = 0;

  /*----- Check for proper number of arguments. -----*/ 
  if(nrhs > 7) { mexErrMsgTxt("up to seven inputs required"); } 
 
  mxGetString(prhs[0], recording_path, NAME_LENGTH);           
  mxGetString(prhs[1], recording_filename_root, NAME_LENGTH);           

  if (nrhs >= 3) {
    decimation_factor = mxGetScalar(prhs[2]); }
  if (nrhs >= 4) {
    comedi_num_channels_to_write = mxGetScalar(prhs[3]);}
  if (nrhs >= 5) {
    nspike_num_channels_to_write_low = mxGetScalar(prhs[4]); }
  if (nrhs >= 5) {
    nspike_num_channels_to_write_high = mxGetScalar(prhs[5]); }
  if (nrhs >= 6) {
    nspike_high_channels_offset = mxGetScalar(prhs[6]); }
    
  NStreamClient *nsc = NStreamClient_create();
    
  if (!NStreamClient_start_record(nsc, (char*)&recording_path, (char*)&recording_filename_root, decimation_factor, comedi_num_channels_to_write, nspike_num_channels_to_write_low, nspike_num_channels_to_write_high, nspike_high_channels_offset))
  {
      NStreamClient_get_error(nsc, (char*)&error_text, NSTREAMCLIENT_ERROR_LENGTH);
      mexErrMsgTxt(error_text);
  }

  NStreamClient_destroy(nsc);
}

/*----------------------------------------------------------------------*/

