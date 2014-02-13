%% Collect parameters for cell model

%{
CharHyperpolStepTA().Baseline_mV
CharHyperpolStepTA().Resistance_Mohms
CharHyperpolStepTA().step_sections
CharHyperpolStepTA().IsThereIhYNeach
CharHyperpolStepTA().output.IhDecayToThird_ms
CharHyperpolStepTA().output.Ih_mV

CharDepolStepTA().Baseline_mV
CharDepolStepTA().Resistance_Mohms
CharDepolStepTA().step_sections
CharDepolStepTA().Spikes_ISIs

CharDepolTonicSpikesTA().ISIs
CharDepolTonicSpikesTA().Spike_ave
CharDepolTonicSpikesTA().SpikeAHP_HalfWidth_ms
CharDepolTonicSpikesTA().Spike_Amplitude
CharDepolTonicSpikesTA().Spike_Width
CharDepolTonicSpikesTA().SpikeRate_persec
%}

funnames = {'CharHyperpolStepTA','CharDepolStepTA','CharDepolTonicSpikesTA'};
%outnames{1} = {'Baseline_mV','Resistance_Mohms',...};

O = CharDepolStepTA(x,y,offset_voltage,tonic_injected_current,sections_label_num,sections_start_sec,sections_length_sec,baseline_start_sec,baseline_length_sec,step_length)
O = CharHyperpolStepTA(x,y,offset_voltage,tonic_injected_current,sections_label_num,sections_start_sec,sections_length_sec,baseline_start_sec,baseline_length_sec)
O = CharDepolTonicSpikesTA(x,yF,yIC,bpFiltParms,Notch,PlotsYN,offset_voltage,tonic_injected_current,varargin)
O = PowerSpecWobbleTA(x,y,FreqRange,Notch,mtmBandWidth,tonic_injected_current,varargin)
O = IPSPsTA(x,yF,yIC,Size,Notch,method,FreqRange,Bins,SmoothWindow,PlotsYN)


%% Network
%{
IPSPsTA().IPSPmeanUni
PowerSpecWobbleTA().????
%}
