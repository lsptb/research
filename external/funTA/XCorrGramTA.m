function XCorrGramTA(x,y1,y2,bpFiltParms,Notch,NormAbs)


% Initial bits n bobs
Fs = fix(1/(x(2)-x(1))); %Fs = round(1/(data{1}(2,1)-data{1}(1,1)));

% Actions
d1 = y1.*1e6;
d2 = y2.*1e6;
if ~isempty(bpFiltParms);
    d1 = Bpfft(d1,Fs,bpFiltParms(1),bpFiltParms(2));
    d2 = Bpfft(d2,Fs,bpFiltParms(1),bpFiltParms(2));
end
if ~isempty(Notch), d1 = Bsfft(d1,Fs,Notch(1),Notch(2)); end
if ~isempty(Notch), d2 = Bsfft(d2,Fs,Notch(1),Notch(2)); end

dim2 = floor(length(d1)/Fs);
d1_mat = reshape(d1(1:dim2*Fs),Fs,dim2);
d2_mat = reshape(d2(1:dim2*Fs),Fs,dim2);
d1_cell = num2cell(d1_mat,1);
d2_cell = num2cell(d2_mat,1);
c = cellfun(@(x,y) xcov(x,y,Fs),d1_cell,d2_cell,'Uni',0);
c2 = c;
c_mat = cell2num(c2);
c_cell = imfilter(interp2(c_mat),fspecial('disk',2));
yomax = max(max((c_cell)));
yomin = min(min((c_cell)));
mmax = max(yomax);
mmin = max(yomin);

if strcmp('Normalized',NormAbs)
    c_cell = (c_cell-mmin)./(mmax-mmin); %c_cell = cellfun(@(x) x./mmax,c_cell,'Uni',0);
end
cmmin = min(min(c_cell));
CMMIN = min(cmmin);

% Calculate axis limits
pk = max(c_cell);
pkmax = max(pk);
mpktemp = max(pkmax);
mpk(2) = mpktemp;
c2 = cellfun(@(x) x*(-1),c,'Uni',0);
pk2 = cellfun(@findpeaks,c2,'Uni',0);
pk2max = cellfun(@max,pk2,'Uni',0);
mpktemp = cellfun(@max,pk2max,'Uni',0);
for k = 1:length(mpktemp), mpktemp2(k) = mpktemp{k}; end
mpktemp = max(mpktemp2);
mpk(1) = mpktemp*(-1);

% Plot results
cmap1 = flipud(cbrewer('div','RdGy',256));
cs = round(size(c_mat,2)/2);
csy = 1000*(-cs/Fs):1000/Fs:1000*(cs/Fs);
fm, imagesc(0.5:0.5:dim2,csy,c_cell');
if strcmp('Normalized',NormAbs)
    caxis([CMMIN 1])
else caxis(mpk)
    ylim([-200 200]), xlabel('Time (secs)'), ylabel('Delay (ms)'), colorbar, colormap(cmap1)
end


end

