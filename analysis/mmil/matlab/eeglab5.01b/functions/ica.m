% ica() - Signal processing functions of the EEGLAB toolbox
%
% TOOLBOX CREDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      EEGLAB -- MATLAB functions for psychophysiological data analysis     %%
%%        including an enhanced and automated version of the infomax         %%
%%             algorithm for Independent Component Analysis (ICA)            %%
%%                        of Bell & Sejnowski (1995).                        %%
%%     by Scott Makeig, Arnaud Delorme, Colin Humphries, Sigurd Enghoff,     %%
%%       with Tzyy-Ping Jung, Tony Bell, Martin McKeown, Luca Finelli        %%
%%          Te-Won Lee, Benjamin Blankertz, Alex Dimitrov, et al.            %%
%%   Swartz Center for Computational Neuroscience, INC, UCSD, Version 4.0    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% GENERAL ELECTROPHYSIOLOGICAL DATA PROCESSING TOOLS:
%
% Find abs peak frames and amplitudes:                            abspeak()
% Change reference from common to average:                        averef()
% Simple block average data epochs:                               blockave()
% Make a 2-D scalp field movie:                                   eegmovie()
% Frequency band filter data:                                     eegfilt()
% View continuous data traces:                                    eegplot()
% Average data epochs (with windowing options):                   erpave()
% Display raw or smoothed single data epochs:         erpimage(), erpimages()
% Re-align event-related epochs to given events:                  eventlock()
% Plot one or more field maps on 3-D head model(s):   headplot(), compheads()
% Construct a movie of a moving field on a 3-D head model:        headmovie()
% Compute and view log power spectra of single data epochs:       logspec()
% Select chans,frames,epochs of concatenated data epochs:         matsel()
% Perform moving averaging on data:                               movav()
% Plot a multichannel data epoch on a single axis:                ploterp()
% Perform principal component analysis (PCA) via SVD              pcasvd()
% Perform nonlinear (post-PCA) rotations:    varimax(), promax(), qrtimax()
% View concatenated multichannel data epochs:                     plotdata()
% View concatenated data epochs in topographic arrangement:       plottopo()
% Plot a data epoch with topoplots at selected time points:       timtopo()
% Change the data sampling rate:                                  resample()
% Remove baseline means from data epochs:                         rmbase()
% View a 2-D or 3-D scalp-field movie:                            seemovie()
% Compute and plot statistics of 1-D data                         signalstat()
% Power spectral scalp topographies                               spectopo()
% Plot time/frequency (ERSP/ITC) scalp topographies               tftopo()
% Event-related time/frequency (ERSP, ITC) of single-trial data:  timef()
% Event-related coherence of single-trial data:                   crossf()
% View data scalp topography(s):                      topoplot(), compmap()
% Convert Cartesian (x,y,z) channel locs to topoplot() format:    cart2topo()
% Convert 2-D topoplot() channel locs to 3-D headplot() format:   topo2sph()
% Convert 2-D headplot() channel locs to 2-D topoplot() format:   sph2topo()
% 
% SPECIFIC ICA TOOLS: 
%
% Perform ICA analysis using logistic infomax or extended-infomax runica() 
% Fastest, most compact: system-call of binary runica()           binica()
% Fast, compact Matlab MEX-file implementation of runica()        mexica() 
% Perform ICA analysis using 2nd & 4th-order cumulants (Cardoso)  jader()
% Test ICA algorithm accuracy, varying data parameters:           testica()
% Compare ICA weight matrices:                       matcorr() -> matperm() 
% Plot data and component envelopes:                    envproj() envtopo() 
% Compute component activations:                                  icaact()
% Compute component variances on scalp:                           icavar()
% Make activations all rms-positive:                              posact()
% Compute component projections:                                  icaproj() 
% Plot the data decomposition:                      plotproj() -> chanproj()
% Plot the data decomposition using plotopo():                    projtopo()
% Sort ICA components by max projected latency and variance:      compsort()
% Sort ICA components by mean projected variance only:            varsort() 
% View a projected ICA component (time course plus topo map):     compplot()
% Squash or expand data into a PCA-defined subspace: pcsquash() -> pcexpand()
% Plot selected time periods of component activations:            tree()
%
% GENERAL HELPER FUNCTIONS:
%
% Make plot axes pop up into zoomable windows on mouse click      axcopy()
% Create possibly-overlapping subplot axes on a general grid      sbplot()
% Plot custom colorbar                                            cbar()
%
% POP-UP TOOLBOX TUTORIAL                                         tutorial()
%
% REFERENCES:          
%    http://sccn.ucsd.edu/eeglab/tutorial/
%    http://sccn.ucsd.edu/eeglab/icabib.html 
% Further information: 
%    http://sccn.ucsd.edu/eeglab/icafaq.html 
%
% SEND news/bugs/fixes/suggestions to: scott@sccn.ucsd.edu

help ica

% $Log: ica.m,v $
% Revision 1.5  2002/11/15 03:20:42  arno
% remove imagetopo function reference
%
% Revision 1.4  2002/11/15 02:34:01  arno
% debug for web
%
% Revision 1.3  2002/08/13 17:01:54  scott
% revised for 4.0
%
% Revision 1.2  2002/04/12 02:59:27  scott
% fixed typo -sm
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% Begun: May 7, 1996
% Current version: Mon Aug 21 13:35:49 PDT 2000
% 01-25-02 reformated help & license, added links -ad 

