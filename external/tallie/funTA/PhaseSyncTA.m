function O = PhaseSyncTA(x,y1,y2,Wind,WindOL)


Fs = fix(1/(x(2)-x(1))); %Fs = round(1/(data{1}(2,1)-data{1}(1,1)));
nx = length(y1);
y1 = y1(:); y2 = y2(:); Wind = Wind(:);
ncol = fix((nx-WindOL)/(Wind-WindOL));
colindex = 1 + (0:(ncol-1))*(Wind-WindOL);
rowindex = (1:Wind)';
z = zeros(Wind+1,ncol);
y = zeros(Wind,ncol);
d = zeros(Wind,ncol);
y(:) = y1(rowindex(:,ones(1,ncol))+colindex(ones(Wind,1),:)-1);
d(:) = y2(rowindex(:,ones(1,ncol))+colindex(ones(Wind,1),:)-1);
next_value = 1;
while next_value < ncol+1;
    q = xcov(y(:,next_value),d(:,next_value),(Wind/2),'coeff');
    z(:,next_value) = q;
    next_value = next_value+1;
end
t = (colindex-1)'/Fs;
fm, imagesc(t,[(-Wind/2)/(Fs/1000) (Wind/2)/(Fs/1000)],(z+eps))
axis xy, caxis([-0.8 0.8]), colormap(jet), xlabel('Time (s)'), ylabel('delay (ms)'), colorbar
sync = z((Wind/2)+1,:);
[~,phase] = max(z,[],1);

O.sync = sync;
O.phase = phase-(Wind/2)+1;
O.raw1 = y1;
O.raw2 = y2;
O.Fs = Fs;
O.msPerPoint = t(2)-t(1);

end

