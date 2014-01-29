function O = SpecGramTA(x,y,Smooth,Notch,Wind,WindOL,SmoothWindow)


Fs = fix(1/(x(2)-x(1))); %Fs = roundn(1/(data{1}(2,1)-data{1}(1,1)),3);
d = y.*1e6;

if ~isempty(Notch), d = Bsfft(d,Fs,Notch(1),Notch(2)); end

switch Smooth
    case {'standard','std'}
        [yo,fo,to] = specmw(detrend(d),Wind,Fs,Wind,WindOL);
    case {'pwelch','pw'}
        [yo,fo,to] = mtmspecTA2(detrend(d),Wind,Fs,Wind,WindOL);
    case {'multitaper','mtm'}
        [yo,fo,to] = mtmspecTA(detrend(d),Wind,Fs,Wind,WindOL,SmoothWindow);
end

fm, surf(to,fo,yo), shading interp, set(gca,'YScale','log','YTick',[[1:10] [15:5:50] [100 200]])
view(2), axis tight, ylim([1 250]), colormap(flipud(colormap('hot'))), colorbar
%fm, imagesc(to,fo,(abs(yo)+eps).^2); axis xy; colormap(flipud(colormap('hot'))), axis tight, ylim([0 250]), colorbar
xlabel('Time (secs)'), ylabel('Frequency (Hz)')

O.raw = y;
O.data = d;
O.method = Smooth;
O.yo = yo;
O.fo = fo;
O.to = to;
O.Fs = Fs;
O.yo_HzPerBin = fo(2)-fo(1);

end

