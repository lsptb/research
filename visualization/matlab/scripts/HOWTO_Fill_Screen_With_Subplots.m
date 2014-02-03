nplots = 8;
nrows  = ceil(sqrt(nplots));
ncols  = nrows;

% nrows  = 8;
% ncols  = 4;

% ------------------------------
dx = 1 / ncols;
dy = 1 / nrows;

xstep = mod((1:ncols)-1,ncols);
ystep = mod((1:nrows)-1,nrows)+1;
xpos  = xstep*dx;
ypos  = 1-ystep*dy;

figure('position',get(0,'screensize'));
for i = 1:nplots
  xi  = mod(i-1,ncols)+1;
  yi  = floor((i-1)./ncols)+1;
  subplot('Position',[xpos(xi) ypos(yi) 1/ncols 1/nrows]); set(gca,'units','normalized');
  text(.5,.5,num2str(i)); axis([0 1 0 1]);
end

