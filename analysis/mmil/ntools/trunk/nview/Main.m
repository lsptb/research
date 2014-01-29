function varargout = Main(varargin)
% MAIN M-file for Main.fig
%      MAIN, by itself, creates a new MAIN or raises the existing
%      singleton*.
%
%      H = MAIN returns the handle to a new MAIN or the handle to
%      the existing singleton*.
%
%      MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN.M with the given input arguments.
%
%      MAIN('Property','Value',...) creates a new MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Main_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Main

% Last Modified by GUIDE v2.5 08-Jan-2010 17:48:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Main_OpeningFcn, ...
                   'gui_OutputFcn',  @Main_OutputFcn, ...
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




% --- Executes just before Main is made visible.
function Main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Main (see VARARGIN)

% Choose default command line output for Main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Main wait for user response (see UIRESUME)
% uiwait(handles.figure1);
%set(0,'Units','pixels') % Sets units to pixels
%set(handles.figure1,'Units','pixels') % Sets units to pixels
%screensize = get(0,'ScreenSize');% Gets exact screen size of your computer
%set(handles.figure1, 'Position',[screensize(1), screensize(2), screensize(3), screensize(4)-75]) % Generates fig. to occupy the whole screen size automatically

% scale to fill screen
units=get(gcf,'units');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'units',units);
set(gcf,'BackingStore', 'off');
set(gcf,'DoubleBuffer', 'on');
set(gcf,'Renderer', 'opengl');

path(path, '/home/nspike/svn_code/nsuite/nstream/pkg/mex');
initialize_global_state()


function initialize_global_state()
% get figure width and height
% units=get(gcf,'units');
% set(gcf,'units','pixels');
% fPosition = get(gcf, 'Position');
% set(gcf,'units',units);


global gAppState experiment

handles = guihandles;

gAppState.recording.state_stopped = 0;
gAppState.recording.state_run     = 1;
gAppState.recording.state = gAppState.recording.state_stopped;
gAppState.recording.recording_path_root = experiment.recording.recording_path;
decimation_factor = round(3e4./experiment.recording.sample_rate);
comedi_decimation_factor = round(3e4./experiment.recording.comedi.sample_rate);
gAppState.recording.decimation_factor = decimation_factor;
gAppState.recording.comedi_decimation_factor = comedi_decimation_factor;
if(3e4 / decimation_factor ~= experiment.recording.sample_rate)
    experiment.recording.sample_rate = 3e4 / decimation_factor;
end
if isfield(experiment.recording,'comedi_num_channels_to_write')
    gAppState.recording.comedi_num_channels_to_write = ...
        experiment.recording.comedi_num_channels_to_write;
else
    gAppState.recording.comedi_num_channels_to_write = 8;
end
if isfield(experiment.recording,'nspike_num_channels_to_write_low')
    gAppState.recording.nspike_num_channels_to_write_low = ...
        experiment.recording.nspike_num_channels_to_write_low;
else
    gAppState.recording.nspike_num_channels_to_write_low = 256;
end
if isfield(experiment.recording,'nspike_num_channels_to_write_high')
	gAppState.recording.nspike_num_channels_to_write_high = ...
        experiment.recording.nspike_num_channels_to_write_high;
else
	gAppState.recording.nspike_num_channels_to_write_high = 0;
end
if isfield(experiment.recording,'nspike_high_channels_offset')
  gAppState.recording.nspike_high_channels_offset = ...
        experiment.recording.nspike_high_channels_offset;
else
  gAppState.recording.nspike_high_channels_offset = 0;
end
  
gAppState.recording.rec_num = 0;

gAppState.traceview.num_traces  = 32;
gAppState.traceview.start_trace = 1;
gAppState.traceview.end_trace   = 32;
gAppState.traceview.vertrange   = 2000;
gAppState.traceview.state       = 0;    % stopped
gAppState.traceview.state_stopped = 0;  % stopped
gAppState.traceview.state_run     = 1;  % running 
gAppState.traceview.state_stop    = 2;  % request stop 
gAppState.traceview.axes          = handles.TraceView_Axes;
gAppState.traceview.sample_rate_num = 3;
gAppState.traceview.sample_rate = 10000;
gAppState.traceview.horiz_time = 1500; % number of milliseconds to display
gAppState.traceview.reference_channel = 0;

gAppState.traceview.filtering_line_noise = 0;
gAppState.traceview.notchfilter = nview_notchfilter([30000 15000 10000 6000 5000 3000 1000]);

for i = 1:256
    gAppState.channels(i).caption = experiment.channels(i).name;
end
gAppState.channels(257).caption = 'Time code';
gAppState.channels(258).caption = 'Audio';
gAppState.channels(259).caption = 'Horiz Eye';
gAppState.channels(260).caption = 'Vert Eye';
gAppState.channels(261).caption = 'Pupil Area';
gAppState.channels(262).caption = 'Photodiode';
gAppState.channels(263).caption = 'Analog 6';
gAppState.channels(264).caption = 'Analog 7';

gAppState.spectralview.running = 0;

build_trace_select_popup(gAppState.traceview.num_traces, handles);
setup_plot_area();

% --- Outputs from this function are returned to the command line.
function varargout = Main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in runbutton.
function runbutton_Callback(hObject, eventdata, handles)
% hObject    handle to runbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global gAppState;
    switch (gAppState.traceview.state)
        case gAppState.traceview.state_stopped 
            gAppState.traceview.state = gAppState.traceview.state_run;
            main_loop();
        case gAppState.traceview.state_run 
            gAppState.traceview.state = gAppState.traceview.state_stop;
    end
   

function main_loop()
    global gAppState;
    setup_plot_area();
    data = zeros(264, 6e4);
    no_ref_data = zeros(1, size(data, 2));
    while(1)
        t = nstream_gettime_nspike;
        t2 = t - gAppState.traceview.horiz_time;
        data = zeros(264,gAppState.traceview.horiz_time/1e3*gAppState.traceview.sample_rate);

        data(1:256,:) = nstream_getrawint_nspike_decimate(t2, t,(30000/gAppState.traceview.sample_rate));
        t_comedi = nstream_gettime_comedi;
        t2_comedi = t_comedi - gAppState.traceview.horiz_time;
        comedi_data = nstream_getrawint_comedi_decimate(t2_comedi, t_comedi,(30000/gAppState.traceview.sample_rate));
        data(257,:) = comedi_data(1,:) - mean(comedi_data(1,:));
        data(258,:) = comedi_data(2,:) - mean(comedi_data(2,:));
        data(259,:) = (comedi_data(3,:)-32768).*.9;
        data(260,:) = (comedi_data(4,:)-32768).*.9;
        data(261,:) = (comedi_data(5,:)-32768).*.9;
        data(262,:) = comedi_data(6,:) - mean(comedi_data(6,:));
        data(263,:) = comedi_data(7,:) - mean(comedi_data(7,:));
        data(264,:) = comedi_data(8,:) - mean(comedi_data(8,:));
        
        if(gAppState.traceview.reference_channel == 0)
            ref_data = no_ref_data(1:size(data, 2));
        else
            ref_data = data(gAppState.traceview.reference_channel,:);
        end
        
        trace_range = gAppState.traceview.start_trace:gAppState.traceview.end_trace;
        
        data(trace_range,:) = bsxfun(@minus, data(trace_range,:), ref_data);
        if(gAppState.traceview.filtering_line_noise)
            data(trace_range,:) = gAppState.traceview.notchfilter(gAppState.traceview.sample_rate_num).filter(data(trace_range,:), 2);
        end

        for i = trace_range
            update_trace(i, data(i,:), gAppState.traceview.traces);
        end

        if(gAppState.spectralview.running)
            spectrumscope(gAppState.spectralview.hax, data(trace_range, :));
        end
        
        if(gAppState.traceview.state == gAppState.traceview.state_stop)
            gAppState.traceview.state = gAppState.traceview.state_stopped;
            break;
        else
            if(gAppState.traceview.filtering_line_noise)
                pause(0.02); % attempt to compensate for increased delay due to filtering
            else
                pause(0.04);
            end
        end
    end
    

function setup_plot_area(height, width, channels, startchan, endchan)
    global gAppState;
    startchan = gAppState.traceview.start_trace %#ok
    endchan   = gAppState.traceview.end_trace %#ok
    channels  = gAppState.channels;
    vertrange = gAppState.traceview.vertrange;
    axes      = gAppState.traceview.axes;
    horzwidth = gAppState.traceview.horiz_time * (gAppState.traceview.sample_rate * .001) %#ok

    plot(zeros(1,horzwidth));
    whitebg('k');
    box on;
    chanheight = vertrange;
    chanoffset = chanheight / 2;
%    axis([0 1e4 -65536 65536]);
%    axis([0 9e3 0 chanheight*(gAppState.traceview.num_traces)]);
    set(gcf, 'Renderer', 'opengl');
    set(axes, 'XLimMode', 'manual');
    set(axes, 'YLimMode', 'manual');
    set(axes, 'XLim', [0 horzwidth]);
    set(axes, 'YLim', [0 chanheight*(gAppState.traceview.num_traces)]);
    set(axes, 'Units', 'pixels');
    %set(axes, 'Position', [50 100 width height]);
    set(axes, 'XTickMode', 'manual');
    set(axes, 'YTickMode', 'manual');
    set(axes, 'XTickLabelMode', 'manual');
    set(axes, 'YTickLabelMode', 'manual');
    set(axes, 'DrawMode', 'fast');
   
    clear ticks; 
    for i = startchan:endchan 
        traces(i).handle = line(1:horzwidth,zeros(1,horzwidth)+chanoffset);
        traces(i).offset = chanoffset;
        ticks(i - (startchan-1)) = chanoffset;
        line(1:horzwidth,zeros(1,horzwidth)+(chanoffset + chanheight/2),'LineStyle',':','Color','white');
        chanoffset = chanoffset+chanheight;
    end

    set(axes, 'YTick', ticks);
    set(axes, 'XTick', []);
    set(axes, 'YTickLabel', {channels(startchan:endchan).caption}); 
    gAppState.traceview.traces = traces;
    gAppState.traceview.chanheight = chanheight;
    
function setup_spectral_view()
    global gAppState;
    
    curr_fig = gcf();
    if(~isfield(gAppState, 'spectralview') || ~isfield(gAppState.spectralview, 'hfig') || ~ishandle(gAppState.spectralview.hfig))
        gAppState.spectralview.hfig = open('Spectral.fig');
        set(gAppState.spectralview.hfig, 'CloseRequestFcn', @close_spectral_view);

        set(findobj(gAppState.spectralview.hfig, 'Tag', 'max_freq'), 'Callback', @update_spectral_freq_window);
        set(findobj(gAppState.spectralview.hfig, 'Tag', 'min_freq'), 'Callback', @update_spectral_freq_window);
    else
        figure(gAppState.spectralview.hfig);
    end

    gAppState.spectralview.hax = gca;
    spectrumscope(gAppState.spectralview.hax, gAppState.traceview.sample_rate, ...
        2^nextpow2(gAppState.traceview.sample_rate * gAppState.traceview.horiz_time / 1000), ...
        gAppState.traceview.num_traces);
    set(gAppState.spectralview.hax, 'Color', 'k', 'XColor', 'y', 'YColor', 'y', 'ZColor', 'y');
    set(findobj(gAppState.spectralview.hax, 'Tag', 'SpectrumScope'), 'Color', 'y');

    hline = findobj(gAppState.spectralview.hax,'Tag','SpectrumScope');
    for i=1:length(hline)
        if(mod(i, 2) == 1)
            line_color = 'yellow';
        else
            line_color = 'cyan';
        end
        set(hline(i), 'Color', line_color);
    end

    update_spectral_freq_window([], []);
    
    figure(curr_fig);
    
function update_spectral_freq_window(src, event)
    global gAppState;

    max_freq = str2double(get(findobj(gAppState.spectralview.hfig, 'Tag', 'max_freq'), 'String'));
    min_freq = str2double(get(findobj(gAppState.spectralview.hfig, 'Tag', 'min_freq'), 'String'));
    if(isnan(max_freq))
        curr_xlim = get(gAppState.spectralview.hax, 'Xlim');
        max_freq = curr_xlim(2);
    end
    if(isnan(min_freq))
        min_freq = 0;
    end
    
    if(min_freq > max_freq)
        disp('Warning: invalid values for spectral frequency limits!');
    end

    xlim(gAppState.spectralview.hax, [min_freq max_freq]);

function close_spectral_view(src, event)
    global gAppState;
    gAppState.spectralview.running = 0;
    delete(gAppState.spectralview.hfig);

function update_trace(trace_number, data, traces)
    line_handle = traces(trace_number).handle;
    set(line_handle, 'YData', data + traces(trace_number).offset);


% --------------------------------------------------------------------
function View_Traces_Show_32_MenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to View_Traces_Show_32_MenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gAppState
gAppState.traceview.start_trace = 1;
gAppState.traceview.end_trace   = 8;
gAppState.traceview.num_traces  = 8;
setup_plot_area();


% --------------------------------------------------------------------
function View_Traces_Show_64_MenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to View_Traces_Show_64_MenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gAppState
gAppState.traceview.start_trace = 1;
gAppState.traceview.end_trace   = 8;
gAppState.traceview.num_traces  = 8;
setup_plot_area();


% --- Executes on selection change in TraceNumber_Menu_PopUp.
function TraceNumber_Menu_PopUp_Callback(hObject, eventdata, handles)
% hObject    handle to TraceNumber_Menu_PopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns TraceNumber_Menu_PopUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TraceNumber_Menu_PopUp
global gAppState
switch get(hObject,'Value')
    % Show 1
    case 1
        gAppState.traceview.start_trace = 1;
        gAppState.traceview.end_trace   = 1;
        gAppState.traceview.num_traces  = 1;
    % Show 2
    case 2
        gAppState.traceview.start_trace = 1;
        gAppState.traceview.end_trace   = 2;
        gAppState.traceview.num_traces  = 2;
    % Show 4
    case 3
        gAppState.traceview.start_trace = 1;
        gAppState.traceview.end_trace   = 4;
        gAppState.traceview.num_traces  = 4;
    % Show 8 
    case 4
        gAppState.traceview.start_trace = 1;
        gAppState.traceview.end_trace   = 8;
        gAppState.traceview.num_traces  = 8;
    % Show 16 
    case 5
        gAppState.traceview.start_trace = 1;
        gAppState.traceview.end_trace   = 16;
        gAppState.traceview.num_traces  = 16;
    % Show 32 
    case 6 
        gAppState.traceview.start_trace = 1;
        gAppState.traceview.end_trace   = 32;
        gAppState.traceview.num_traces  = 32;
    otherwise
end
setup_plot_area();
set(handles.Trace_Select,'Value',1);
build_trace_select_popup(gAppState.traceview.num_traces, handles);
if(gAppState.spectralview.running)
    setup_spectral_view();
end

function build_trace_select_popup(increment, handles)
global gAppState
val = {};
for i = 1:increment:264
    if(increment > 1)
        val{end+1} = [gAppState.channels(i).caption '-' gAppState.channels(min(end,i+increment-1)).caption];
    else
        val{end+1} = gAppState.channels(i).caption;
    end
end
   
set(handles.Trace_Select, 'String', val);  
    
% --- Executes during object creation, after setting all properties.
function TraceNumber_Menu_PopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TraceNumber_Menu_PopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function TraceNumber_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to TraceNumber_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('tracenumber menu callback');

% --------------------------------------------------------------------
function Show_2_Traces_MenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to Show_4_Traces_MenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gAppState
gAppState.traceview.start_trace = 1;
gAppState.traceview.end_trace   = 2;
gAppState.traceview.num_traces  = 2;
setup_plot_area();

% --------------------------------------------------------------------
function Show_4_Traces_MenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to Show_4_Traces_MenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gAppState
gAppState.traceview.start_trace = 1;
gAppState.traceview.end_trace   = 4;
gAppState.traceview.num_traces  = 4;
setup_plot_area();


% --------------------------------------------------------------------
function Show_8_Traces_MenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to Show_8_Traces_MenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gAppState
gAppState.traceview.start_trace = 1;
gAppState.traceview.end_trace   = 8;
gAppState.traceview.num_traces  = 8;
setup_plot_area();


% --------------------------------------------------------------------
function Show_16_Traces_MenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to Show_16_Traces_MenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gAppState
gAppState.traceview.start_trace = 1;
gAppState.traceview.end_trace   = 16;
gAppState.traceview.num_traces  = 16;
setup_plot_area();


% --------------------------------------------------------------------
function Show_32_Traces_MenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to Show_32_Traces_MenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gAppState
gAppState.traceview.start_trace = 1;
gAppState.traceview.end_trace   = 32;
gAppState.traceview.num_traces  = 32;
setup_plot_area();
disp('32');


% --- Executes on selection change in Trace_Select.
function Trace_Select_Callback(hObject, eventdata, handles)
% hObject    handle to Trace_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Trace_Select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Trace_Select
global gAppState
selection = get(hObject,'Value');
gAppState.traceview.start_trace = (selection - 1)*gAppState.traceview.num_traces + 1;
gAppState.traceview.end_trace = min(264, gAppState.traceview.start_trace + (gAppState.traceview.num_traces - 1));
gAppState.traceview.start_trace
gAppState.traceview.end_trace
setup_plot_area();

% --- Executes during object creation, after setting all properties.
function Trace_Select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Trace_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Vertical_Range_PopUpMenu.
function Vertical_Range_PopUpMenu_Callback(hObject, eventdata, handles)
% hObject    handle to Vertical_Range_PopUpMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Vertical_Range_PopUpMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Vertical_Range_PopUpMenu
global gAppState
switch get(hObject,'Value')
    % +/- 5000 
    case 1
        gAppState.traceview.vertrange = (5000/10)*2;
    % +/- 10000 
    case 2
        gAppState.traceview.vertrange = (10000/10)*2;
    % +/- 15000 
    case 3
        gAppState.traceview.vertrange = (15000/10)*2;
    % +/- 30000 
    case 4 
        gAppState.traceview.vertrange = (30000/10)*2;
    % +/- 60000 
    case 5 
        gAppState.traceview.vertrange = (60000/10)*2;
    %  Comedi
    case 6
        gAppState.traceview.vertrange = 32768*2;
    otherwise
end
setup_plot_area();


% --- Executes during object creation, after setting all properties.
function Vertical_Range_PopUpMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vertical_Range_PopUpMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton50.
function pushbutton50_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton51.
function pushbutton51_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in spectral_view.
function spectral_view_Callback(hObject, eventdata, handles)
% hObject    handle to spectral_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gAppState;
if(gAppState.spectralview.running)
    gAppState.spectralview.running = 0;
    close_spectral_view([], []);
else
    gAppState.spectralview.running = 1;
    setup_spectral_view();
    figure(gAppState.spectralview.hfig);
end


% --- Executes on button press in pushbutton53.
function pushbutton53_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton54.
function pushbutton54_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton55.
function pushbutton55_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton56.
function pushbutton56_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton57.
function pushbutton57_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Recording_PushButton.
function Recording_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Recording_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gAppState;
switch gAppState.recording.state
    case gAppState.recording.state_stopped
        recordingfilepath = start_recording();
        if(recordingfilepath ~= 0)
            gAppState.recording.state = gAppState.recording.state_run;
            set(hObject, 'String', 'Stop Recording');
            set(handles.stRecordingFileLocation,'String',['Recording data to: ' recordingfilepath]);
            set(handles.RecordingNumberTextArea, 'String', gAppState.recording.rec_num);
        end
    case gAppState.recording.state_run
        set(hObject, 'String', 'Start Recording');
        gAppState.recording.state = gAppState.recording.state_stopped;
        stop_recording();
        set(handles.stRecordingFileLocation,'String',['Stopped recording data']);
end

function rec_num = next_rec_num(path_root)
    old_dir = pwd();
    cd(path_root);
    rec_num_dirs = filenames('[0-9][0-9][0-9]');
    if(isempty(rec_num_dirs))
        rec_num = 1;
    else
        rec_num = str2double(rec_num_dirs{end}) + 1;
    end
    cd(old_dir)

function recording_path = start_recording()
    global gAppState experiment;
    
    gAppState.recording.rec_num = next_rec_num(gAppState.recording.recording_path_root);
    gAppState.recording.rec_dir = sprintf('%03d', gAppState.recording.rec_num);
    new_path = [gAppState.recording.recording_path_root '/' gAppState.recording.rec_dir];
    
    mkdir(new_path);
    gAppState.recording.rec_fileroot = ['rec' sprintf('%03d', gAppState.recording.rec_num)];
    [gAppState.recording.rec_fileroot, recording_path] = ...
    	uiputfile('', 'Save recording as...', fullfile(new_path, gAppState.recording.rec_fileroot));
    	
    if(recording_path ~= 0)
        overwrite = 'No';
        while(strcmp(overwrite, 'No') && exist(fullfile(recording_path, [gAppState.recording.rec_fileroot '.nspike.dat']), 'file'))
        	overwrite = questdlg(sprintf('There appear to already be files named after %s. Are you sure you want to continue?', gAppState.recording.rec_fileroot), 'Confirm overwrite', 'Yes', 'No', 'No');
            if(strcmp(overwrite, 'No'))
                [gAppState.recording.rec_fileroot, recording_path] = ...
                    uiputfile('', 'Save recording as...', fullfile(recording_path, gAppState.recording.rec_fileroot));
            end
        end
        
        % Copy initialization files for archival and start diary
        diary(fullfile(recording_path, 'matlab.log'));
        copyfile(experiment.initialization.patient_defn_file, recording_path);
        copyfile(experiment.initialization.hw_defn_file, recording_path);
        copyfile(experiment.initialization.exp_defn_file, recording_path);
        copyfile(experiment.initialization.montage_file, recording_path);
    
        recording = gAppState.recording; %just to cut down on visual noise
        nstream_start_record(recording_path, recording.rec_fileroot, recording.decimation_factor, ...
            recording.comedi_num_channels_to_write, recording.nspike_num_channels_to_write_low, ...
            recording.nspike_num_channels_to_write_high, recording.nspike_high_channels_offset);
        
    else
        rmdir(new_path);
    end

function stop_recording()
    global gAppState experiment; %#ok;
    MAX_SAMPLING_FREQ = 30000;

    nstream_stop_record();
    
    recording = gAppState.recording;
    recording_path = fullfile(recording.recording_path_root, recording.rec_dir);
        
    full_file_root = fullfile(recording_path, gAppState.recording.rec_fileroot);
    mat_file = [full_file_root '-exp.mat'];
    disp(['Saving experiment to: ' mat_file]);
    save(mat_file, 'experiment');

    % Loop, checking file size until nstream has flushed out to disk
    low_dat_filename = fullfile(recording_path, [recording.rec_fileroot '.low.nspike.dat']);
    file_sizes = [file_size(low_dat_filename) zeros(1, 4)];
    tic;
    while(length(unique(file_sizes)) > 1 && toc < 60)
        pause(2);
        file_sizes = [file_size(low_dat_filename) file_sizes(1:4)];
    end

    if(recording.nspike_num_channels_to_write_high > 0)
        high_channel_range = 1:recording.nspike_num_channels_to_write_high + recording.nspike_high_channels_offset;
    else
        high_channel_range = [];
    end
    ntools_hdrxml_write(fullfile(recording_path, [recording.rec_fileroot '.low.nspike.hdrxml']), MAX_SAMPLING_FREQ / recording.decimation_factor, recording.nspike_num_channels_to_write_low, ...
        [], 'int16', {gAppState.channels(1:recording.nspike_num_channels_to_write_low).caption});
    
    ntools_hdrxml_write(fullfile(recording_path, [recording.rec_fileroot '.high.nspike.hdrxml']), MAX_SAMPLING_FREQ, recording.nspike_num_channels_to_write_high, ...
        [], 'int16', {gAppState.channels(high_channel_range).caption});

    ntools_hdrxml_write(fullfile(recording_path, [recording.rec_fileroot '.comedi.hdrxml']), MAX_SAMPLING_FREQ / recording.comedi_decimation_factor, recording.comedi_num_channels_to_write, ...
        [], 'int16', {gAppState.channels(257:257+recording.comedi_num_channels_to_write-1).caption});
    
    msgbox('Successfully finished writing files to disk.')
    




% --- Executes on button press in pushbutton59.
function pushbutton59_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in sample_rate.
function sample_rate_Callback(hObject, eventdata, handles)
% hObject    handle to sample_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns sample_rate contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sample_rate
global gAppState
gAppState.traceview.sample_rate_num = get(hObject,'Value');
switch gAppState.traceview.sample_rate_num
    % 30khz 
    case 1
        gAppState.traceview.sample_rate = 30000;
    % 15khz 
    case 2
        gAppState.traceview.sample_rate = 15000;
    % 10khz 
    case 3
        gAppState.traceview.sample_rate = 10000;
    % 6khz 
    case 4 
        gAppState.traceview.sample_rate = 6000;
    % 5khz 
    case 5 
        gAppState.traceview.sample_rate = 5000;
    % 3khz 
    case 6 
        gAppState.traceview.sample_rate = 3000;
    % 1khz 
    case 7 
        gAppState.traceview.sample_rate = 1000;
    otherwise
end
setup_plot_area();
gAppState.traceview.sample_rate
if(gAppState.spectralview.running)
    setup_spectral_view();
end


% --- Executes during object creation, after setting all properties.
function sample_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sample_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
%nstream_exit;
delete(hObject);
%pause(1.5);
%exit;



function DisplayWidthTextArea_Callback(hObject, eventdata, handles)
% hObject    handle to DisplayWidthTextArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DisplayWidthTextArea as text
%        str2double(get(hObject,'String')) returns contents of DisplayWidthTextArea as a double
global gAppState;
input = str2double(get(hObject,'String'));
 
%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String', gAppState.traceview.horiz_time)
     guidata(hObject, handles);
     return;
end

gAppState.traceview.horiz_time = input;
setup_plot_area();
if(gAppState.spectralview.running)
    setup_spectral_view();
end

% --- Executes during object creation, after setting all properties.
function DisplayWidthTextArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DisplayWidthTextArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on DisplayWidthTextArea and none of its controls.
function DisplayWidthTextArea_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to DisplayWidthTextArea (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)



function edDAC_Channel_0_Callback(hObject, eventdata, handles)
% hObject    handle to edDAC_Channel_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edDAC_Channel_0 as text
%        str2double(get(hObject,'String')) returns contents of edDAC_Channel_0 as a double

Channel = str2num(get(handles.edDAC_Channel_0,'String'));
nstream_set_dac_channel(0,Channel);

% --- Executes during object creation, after setting all properties.
function edDAC_Channel_0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edDAC_Channel_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edDAC_Channel_1_Callback(hObject, eventdata, handles)
% hObject    handle to edDAC_Channel_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edDAC_Channel_1 as text
%        str2double(get(hObject,'String')) returns contents of edDAC_Channel_1 as a double

Channel = str2num(get(handles.edDAC_Channel_1,'String'));
nstream_set_dac_channel(1,Channel);



% --- Executes during object creation, after setting all properties.
function edDAC_Channel_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edDAC_Channel_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edDAC_Gain_0_Callback(hObject, eventdata, handles)
% hObject    handle to edDAC_Gain_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edDAC_Gain_0 as text
%        str2double(get(hObject,'String')) returns contents of edDAC_Gain_0 as a double

Gain = str2num(get(handles.edDAC_Gain_0,'String'));
nstream_set_dac_gain(0,Gain);

% --- Executes during object creation, after setting all properties.
function edDAC_Gain_0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edDAC_Gain_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edDAC_Gain_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edDAC_Gain_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edDAC_Gain_1_Callback(hObject, eventdata, handles)
% hObject    handle to edDAC_Gain_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edDAC_Gain_1 as text
%        str2double(get(hObject,'String')) returns contents of edDAC_Gain_1 as a double

Gain = str2num(get(handles.edDAC_Gain_1,'String'));
nstream_set_dac_gain(1,Gain);


% --- Executes on selection change in RefChannel_Popup_Menu.
function RefChannel_Popup_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to RefChannel_Popup_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns RefChannel_Popup_Menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RefChannel_Popup_Menu
global gAppState

gAppState.traceview.reference_channel = get(hObject,'Value') - 1;
%if(~strcmp(channel_selection, 'None'))
%    channel_selection = str2double(channel_selection);
%    gAppState.traceview.reference_channel = channel_selection;
%else
%    gAppState.traceview.reference_channel = 0;
%end

% --- Executes during object creation, after setting all properties.
function RefChannel_Popup_Menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RefChannel_Popup_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

global experiment;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%for i=1:length(experiment.channels)
%    experiment.channels(i).name
%end

set(hObject, 'String', {'None', experiment.channels(:).name});


% --- Executes during object creation, after setting all properties.
function TraceView_Axes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TraceView_Axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate TraceView_Axes


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in line_noise_filter.
function line_noise_filter_Callback(hObject, eventdata, handles)
% hObject    handle to line_noise_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of line_noise_filter
global gAppState

gAppState.traceview.filtering_line_noise = (get(hObject,'Value'));


