function TA_SlideDisplay_NEW(PlotType,x_name,y_name,Z_name,varargin)

% N.B. No figure is created.  If a figure already exists, this will be
% overlaided on top of it.
%
% PlotType  = either 'plot' or 'imagesc'.
% x_name    = name of the xdata variable in the base workspace to plot.
% y_name    = name of the ydata variable in the base workspace to plot.
% Z_name    = name of the 2D matrix in the base workspace to image - LEAVE
%             AS AN EMPTY MATRIX IF NOT PLOTTING AN IMAGE.
% varargin paired inputs, i.e. Property Name and Property Value:
% (up to three paired inputs are allowed)
%   SlideFun = cell array of strings specifying the function handle for the
%              sliders and it's inputs.  The input you want to remain
%              variable (i.e. the slider value) should be an empty matrix.
%              If using multiple sliders, list the input the functions
%              in the order that you want them to happen to the data.
%   SlideLim = A two-element vector of min and max values for the sliders.
%
% Functions tested (or created for this) so far:
% interp2
% imfilterSubFun
% downsample
% sgolayfilt
% Bpfft

% create sliders
var_num = numel(varargin);
    SlideFunCell{1} = ['1' varargin{1}]; SlideFun{1} = SlideFunCell{1}{2}; SlideLim = varargin{2};
    hSlid1 = uicontrol('Style','Slider','Min',SlideLim(1),'Max',SlideLim(2),'Value',round(SlideLim(2)/4),'Position',[10 10 400 15],'Callback',{@SUB1});
    uicontrol('Style','text','Position',[410 10 200 15],'String',[SlideFun{1} ' ' num2str(get(hSlid1,'Value'))]);
if var_num>2
    SlideFunCell{2} = ['1' varargin{3}]; SlideFun{2} = SlideFunCell{2}{2}; SlideLim = varargin{4};
    hSlid2 = uicontrol('Style','Slider','Min',SlideLim(1),'Max',SlideLim(2),'Value',round(SlideLim(2)/4),'Position',[10 30 400 15],'Callback',{@SUB2});
    uicontrol('Style','text','Position',[410 30 200 15],'String',[SlideFun{2} ' ' num2str(get(hSlid2,'Value'))]);
end
if var_num>4
    SlideFunCell{3} = ['1' varargin{5}]; SlideFun{3} = SlideFunCell{3}{2}; SlideLim = varargin{6};
    hSlid3 = uicontrol('Style','Slider','Min',SlideLim(1),'Max',SlideLim(2),'Value',round(SlideLim(2)/4),'Position',[10 50 400 15],'Callback',{@SUB3});
    uicontrol('Style','text','Position',[410 50 200 15],'String',[SlideFun{3} ' ' num2str(get(hSlid3,'Value'))]);
end

% import workspace vars and perform functions with default slider settings
input(1) = get(hSlid1,'Value');
if exist('hSlid2'), input(2) = get(hSlid2,'Value'); else hSlid2 = []; end
if exist('hSlid3'), input(3) = get(hSlid3,'Value'); else hSlid3 = []; end
switch PlotType
    case 'plot'
        OUT.out1 = evalin('base',[y_name ';']); assignin('base','OUT',OUT);
        for c1 = 1:length(SlideFun)
            [~,j1(c1)] = find(cellfun(@isempty,SlideFunCell{c1}));
            SlideFunCell{c1}(j1(c1)) = {num2str((round(input(c1))))};
            [~,j2(c1)] = find(cellfun(@(x) strcmp(y_name,x),SlideFunCell{c1}));
            SlideFunCell{c1}(j2(c1)) = {'OUT.out1'};
            OUT = EvalCellFun(SlideFunCell{c1}); assignin('base','OUT',OUT);
        end
    case 'imagesc'
        OUT.out1 = evalin('base',[Z_name ';']); assignin('base','OUT',OUT);
        for c1 = 1:length(SlideFun)
            [~,j1(c1)] = find(cellfun(@isempty,SlideFunCell{c1}));
            SlideFunCell{c1}(j1(c1)) = {num2str((round(input(c1))))};
            [~,j2(c1)] = find(cellfun(@(x) strcmp(Z_name,x),SlideFunCell{c1}));
            SlideFunCell{c1}(j2(c1)) = {'OUT.out1'};
            OUT = EvalCellFun(SlideFunCell{c1}); assignin('base','OUT',OUT);
        end
end

switch PlotType
    case 'plot'
        SUBplot(x_name,y_name,SlideFun,input,1,hSlid1,hSlid2,hSlid3)
        ym = evalin('base',['nanmedian(' y_name ');']); om = nanmedian(OUT.out1); ymom = ym-om;
        hold on, hp = plot(xt,OUT.out1+ymom,'-b');
        uicontrol('Style','text','Position',[410 10 200 15],'String',[SlideFun{1} ' ' num2str(get(hSlid1,'Value'))]);
    case 'imagesc'
        SUBimagesc(x_name,y_name,Z_name,SlideFun,input,1,hSlid1,hSlid2,hSlid3)
        hold on, hp = imagesc(xt,yt,OUT.out1);
        uicontrol('Style','text','Position',[410 10 200 15],'String',[SlideFun{1} ' ' num2str(get(hSlid1,'Value'))]);
        axis tight, axis xy
end



% NESTED FUNCTIONS

    function SUB1(hSlid1,event)
        input(1) = get(hSlid1,'Value');
        switch PlotType
            case 'plot'
                SUBplot(x_name,y_name,SlideFun,input,1,hSlid1,hSlid2,hSlid3)
                ym = evalin('base',['nanmedian(' y_name ');']); om = nanmedian(OUT.out1); ymom = ym-om;
                hold on, hp = plot(xt,OUT.out1+ymom,'-b');
                uicontrol('Style','text','Position',[410 10 200 15],'String',[SlideFun{1} ' ' num2str(get(hSlid1,'Value'))]);
            case 'imagesc'
                SUBimagesc(x_name,y_name,Z_name,SlideFun,input,1,hSlid1,hSlid2,hSlid3)
                hold on, hp = imagesc(xt,yt,OUT.out1);
                uicontrol('Style','text','Position',[410 10 200 15],'String',[SlideFun{1} ' ' num2str(get(hSlid1,'Value'))]);
        end
    end
    function SUB2(hSlid2,event)
        input(2) = get(hSlid2,'Value');
        switch PlotType
            case 'plot'
                SUBplot(x_name,y_name,SlideFun,input,1,hSlid1,hSlid2,hSlid3)
                ym = evalin('base',['nanmedian(' y_name ');']); om = nanmedian(OUT.out1); ymom = ym-om;
                hold on, hp = plot(xt,OUT.out1+ymom,'-b');
                uicontrol('Style','text','Position',[410 30 200 15],'String',[SlideFun{2} ' ' num2str(get(hSlid2,'Value'))]);
            case 'imagesc'
                SUBimagesc(x_name,y_name,Z_name,SlideFun,input,1,hSlid1,hSlid2,hSlid3)
                hold on, hp = imagesc(xt,yt,OUT.out1);
                uicontrol('Style','text','Position',[410 30 200 15],'String',[SlideFun{2} ' ' num2str(get(hSlid2,'Value'))]);
        end
    end
    function SUB3(hSlid3,event)
        input(3) = get(hSlid3,'Value');
        switch PlotType
            case 'plot'
                SUBplot(x_name,y_name,SlideFun,input,1,hSlid1,hSlid2,hSlid3)
                ym = evalin('base',['nanmedian(' y_name ');']); om = nanmedian(OUT.out1); ymom = ym-om;
                hold on, hp = plot(xt,OUT.out1+ymom,'-b');
                uicontrol('Style','text','Position',[410 50 200 15],'String',[SlideFun{3} ' ' num2str(get(hSlid3,'Value'))]);
            case 'imagesc'
                SUBimagesc(x_name,y_name,Z_name,SlideFun,input,1,hSlid1,hSlid2,hSlid3)
                hold on, hp = imagesc(xt,yt,OUT.out1);
                uicontrol('Style','text','Position',[410 50 200 15],'String',[SlideFun{3} ' ' num2str(get(hSlid3,'Value'))]);
        end
    end

    function SUBplot(x_name,y_name,SlideFun,input,Slidn,hSlid1,hSlid2,hSlid3)
        OUT.out1 = evalin('base',[y_name ';']); assignin('base','OUT',OUT);
        xt = evalin('base',[x_name ';']);
        for cn = 1:length(SlideFun)
            SlideFunCell{cn}(j1(cn)) = {num2str((round(input(cn))))};
            SlideFunCell{cn}(j2(cn)) = {'OUT.out1'};
            OUT = EvalCellFun(SlideFunCell{cn});
            assignin('base','OUT',OUT);
        end
        if exist('hp'), delete(hp), end %#ok<*EXIST>
        if length(xt)-length(OUT.out1)==1
            xt(1) = [];
        elseif length(OUT.out1)-length(xt)==1
            OUT.out1(1) = [];
        elseif any(strcmp('downsample',SlideFun)) || any(strcmp('decimate',SlideFun));
            d = find(strcmp('downsample',SlideFun));
            xt = downsample(xt,round(input(d))); %#ok<*FNDSB>
        elseif ~mod(length(OUT.out1)/length(xt),1)
            xt = linspace(xt(1),xt(end),length(OUT.out1));
        end
        assignin('base','xt',xt);
    end
    function SUBimagesc(x_name,y_name,Z_name,SlideFun,input,Slidn,hSlid1,hSlid2,hSlid3) %#ok<*INUSD>
        OUT.out1 = evalin('base',[Z_name ';']); assignin('base','OUT',OUT);
        xt = evalin('base',[x_name ';']);
        yt = evalin('base',[y_name ';']);
        for cn = 1:length(SlideFun)
            SlideFunCell{cn}(j1(cn)) = {num2str((round(input(cn))))};
            SlideFunCell{cn}(j2(cn)) = {'OUT.out1'};
            OUT = EvalCellFun(SlideFunCell{cn}); assignin('base','OUT',OUT);
        end
        if exist('hp'), delete(hp), end
        if length(xt)-size(OUT.out1,1)==1
            xt(1) = [];
        elseif length(xt)-size(OUT.out1,1)==(-1)
            OUT.out1(1) = [];
        elseif any(strcmp('downsample',SlideFun)) || any(strcmp('decimate',SlideFun));
            d = find(strcmp('downsample',SlideFun));
            xt = downsample(xt,round(input(d)));
        elseif ~mod(size(OUT.out1,1)/length(xt),1)
            xt = linspace(xt(1),xt(end),length(OUT.out1));
        end
        if length(yt)-size(OUT.out1,1)==1
            yt(1) = [];
        elseif length(yt)-size(OUT.out1,1)==(-1)
            OUT.out1(1) = [];
        elseif any(strcmp('downsample',SlideFun)) || any(strcmp('decimate',SlideFun));
            d = find(strcmp('downsample',SlideFun));
            yt = downsample(yt,round(input(d)));
        elseif ~mod(size(OUT.out1,1)/length(yt),1)
            yt = linspace(yt(1),yt(end),length(OUT.out1));
        end
    end


end

