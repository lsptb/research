function varargout = NTools_Process(varargin)
% NTOOLS_PROCESS M-file for NTools_Process.fig
%      NTOOLS_PROCESS, by itself, creates a new NTOOLS_PROCESS or raises the existing
%      singleton*.
%
%      H = NTOOLS_PROCESS returns the handle to a new NTOOLS_PROCESS or the handle to
%      the existing singleton*.
%
%      NTOOLS_PROCESS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NTOOLS_PROCESS.M with the given input arguments.
%
%      NTOOLS_PROCESS('Property','Value',...) creates a new NTOOLS_PROCESS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NTools_Process_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NTools_Process_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NTools_Process

% Last Modified by GUIDE v2.5 06-May-2009 11:46:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NTools_Process_OpeningFcn, ...
                   'gui_OutputFcn',  @NTools_Process_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before NTools_Process is made visible.
function NTools_Process_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NTools_Process (see VARARGIN)

% Choose default command line output for NTools_Process
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NTools_Process wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NTools_Process_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in convert.
function convert_Callback(hObject, eventdata, handles)
% hObject    handle to convert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [filename, pathname] = uigetfile('*.dat', 'Select a .dat file to convert'); if(filename == 0), return; end
    
    sfreqh = findobj(gcf, 'Tag', 'sfreq');
    sfreq_list = get(sfreqh, 'String');
    new_sfreq = str2double(sfreq_list{get(sfreqh, 'Value')});
    
    precisionh = findobj(gcf, 'Tag', 'precision');
    precision_list = get(precisionh, 'String');
    new_precision = precision_list{get(precisionh, 'Value')};
    
    ntools_proc_convert(fullfile(pathname, filename), new_precision, new_sfreq);
    
    [filename, pathname] = uigetfile('*.dio.txt', 'Select a .dio.txt file to convert or hit cancel to skip.'); if(filename == 0), return; end
    ntools_dio_downsample(fullfile(pathname, filename), new_sfreq);


% --- Executes on button press in filter_line_freqs.
function filter_line_freqs_Callback(hObject, eventdata, handles)
% hObject    handle to filter_line_freqs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global experiment
if(isempty(experiment))
    warnh = warndlg('experiment variable is not set. Please press button 1 to select a .mat file containing the experiment variable.', 'ExperimentWarning', 'modal');
    uiwait(warnh);
else
    [filename, pathname] = uigetfile('*.ieeg.dat', 'Select a .ieeg.dat file');
    ntools_procCleanIEEG(fullfile(pathname, filename));
end

% --- Executes on button press in gen_epoch_data.
function gen_epoch_data_Callback(hObject, eventdata, handles)
% hObject    handle to gen_epoch_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [dat_filename, dat_pathname] = uigetfile('*.dat', 'Select a .dat file'); if(dat_filename == 0), return; end
    [dio_filename, dio_pathname] = uigetfile('*.dio.txt', 'Select a .dio.txt file containing the DIO info'); if(dio_filename == 0), return; end
    [exp_filename, exp_pathname] = uigetfile('*.mat', 'Select a .mat file containing the experiment variable'); if(exp_filename == 0), return; end
    load(fullfile(exp_pathname, exp_filename));
    epoch_data = ntools_gen_epoch_data(fullfile(dat_pathname, dat_filename), fullfile(dio_pathname, dio_filename), experiment, 'coords', []); %#ok

    [d recording_filename_root] = fileparts(dat_filename); %#ok
    [filename, pathname] = uiputfile('*.mat', 'Save file as', [recording_filename_root '.epoch_data.mat']);
    save(fullfile(pathname, filename), '-v7.3', 'epoch_data');


% --- Executes on selection change in sfreq.
function sfreq_Callback(hObject, eventdata, handles)
% hObject    handle to sfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns sfreq contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sfreq


% --- Executes during object creation, after setting all properties.
function sfreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in display_info.
function display_info_Callback(hObject, eventdata, handles)
% hObject    handle to display_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [filename, pathname] = uigetfile({'*.dat';'*.hdrxml'}, 'Select a .dat or .hdrxml file');  if(filename == 0), return; end
    [sfreq, num_channels, num_samples, precision, chan_names] = ntools_hdrxml_read(fullfile(pathname, filename));
    set(findobj(gcf, 'Tag', 'filename_text'), 'String', filename);
    set(findobj(gcf, 'Tag', 'sfreq_text'), 'String', sfreq);
    set(findobj(gcf, 'Tag', 'duration_text'), 'String', [num2str(num_samples/sfreq) ' s']);
    set(findobj(gcf, 'Tag', 'num_channels_text'), 'String', num_channels);
    set(findobj(gcf, 'Tag', 'data_type_text'), 'String', precision);


% --- Executes on selection change in precision.
function precision_Callback(hObject, eventdata, handles)
% hObject    handle to precision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns precision contents as cell array
%        contents{get(hObject,'Value')} returns selected item from precision


% --- Executes during object creation, after setting all properties.
function precision_CreateFcn(hObject, eventdata, handles)
% hObject    handle to precision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gen_ft_data.
function gen_ft_data_Callback(hObject, eventdata, handles)
% hObject    handle to gen_ft_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [dat_filename, dat_pathname] = uigetfile('*.dat', 'Select a .dat file'); if(dat_filename == 0), return; end
    [dio_filename, dio_pathname] = uigetfile('*.dio.txt', 'Select a .dio.txt file containing the DIO info'); if(dio_filename == 0), return; end
    [exp_filename, exp_pathname] = uigetfile('*.mat', 'Select a .mat file containing the experiment variable'); if(exp_filename == 0), return; end
    load(fullfile(exp_pathname, exp_filename));
    ft_data = ntools_gen_ft_data(fullfile(dat_pathname, dat_filename), fullfile(dio_pathname, dio_filename), experiment); %#ok

    [d recording_filename_root] = fileparts(dat_filename); %#ok
    [filename, pathname] = uiputfile('*.mat', 'Save file as', [recording_filename_root '.ft_data.mat']);
    save(fullfile(pathname, filename), '-v7.3', 'ft_data');


% --- Executes during object creation, after setting all properties.
function gen_ft_data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gen_ft_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


