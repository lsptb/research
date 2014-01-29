%     SOreference(refcnt).reflabel
%     SOreference(refcnt).refpeaktype
%     SOreference(refcnt).clusters
%                                 (pktype).CohPSD_corrcoef
%                                 (pktype).CohPSD_corrcoef_p
%     SOreference(refcnt).events
%     SOreference(refcnt).peaks
%     SOreference(refcnt).flip_vec

sel1  = ~arrayfun(@(x)(isempty(x.clusters(1).CohPSD_corrcoef)),SOreference);
sel2  = ~arrayfun(@(x)(isempty(x.clusters(2).CohPSD_corrcoef)),SOreference);

Rneg  = arrayfun(@(x)(x.clusters(1).CohPSD_corrcoef),SOreference(sel1));
Rpos  = arrayfun(@(x)(x.clusters(2).CohPSD_corrcoef),SOreference(sel2));
pneg  = arrayfun(@(x)(x.clusters(1).CohPSD_corrcoef_p),SOreference(sel1));
ppos  = arrayfun(@(x)(x.clusters(2).CohPSD_corrcoef_p),SOreference(sel2));

R     = cellfun(@(x,y)[x y],num2cell(Rneg),num2cell(Rpos),'UniformOutput',false);
p     = cellfun(@(x,y)[x y],num2cell(pneg),num2cell(ppos),'UniformOutput',false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
neglop = cellfun(@(x)x(1)<x(2),p); %% compare to rows of the flip matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rsort = cellfun(@(x)(sort(x)),R,'UniformOutput',false);
psort = cellfun(@(x)(sort(x)),p,'UniformOutput',false);

Rlo   = cellfun(@(x)x(1),Rsort);
Rhi   = cellfun(@(x)x(2),Rsort);
plo   = cellfun(@(x)x(1),psort);
phi   = cellfun(@(x)x(2),psort);

ix1   = 1:length(R);
ix2   = 1:length(p);

figure('Name','Mean PSD vs Coherence (30-50Hz)');
subplot(1,2,1),plot(ix1,Rlo,'r.-',ix1,Rhi,'b.-');legend('lo','hi');title('R(coh,psd)');axis tight
subplot(1,2,2),plot(ix2,plo,'r.-',ix2,phi,'b.-');legend('lo','hi');title('p-value');axis tight;hline(.05,'k');hline(.25,'k');

pdiv = cellfun(@(x)x(2)/x(1),psort);
figure,plot(pdiv,'.');title('maxp/minp'); axis tight; hline(3,'k');

allp = [p{:}];
allR = [R{:}];
figure('Name','Mean PSD vs Coherence for both pos and neg peaks');
subplot(1,3,1),try hist(allR,max(3,floor(length(allR)/10))); end; title('R histogram')
subplot(1,3,2),try hist(allp,max(3,floor(length(allR)/10))); end; title('p-value histogram')
subplot(1,3,3),try hist(pdiv,max(3,floor(length(allR)/10))); end; title('pdiv histogram')

%%

cc=1; 
for k=1:size(flip(cc).matrix,1)
  flipvec = flip(cc).matrix(k,sel2); 
  flipvec(flipvec==-1)=0; 
  rr=corrcoef(neglop,flipvec); 
  r(k)=rr(1,2); 
end; 
figure; 
subplot(1,2,1),plot(r,'.'),axis tight, 
subplot(1,2,2),hist(r,60)

thresh = .03;
rg1 = r(1:2:204);
rg2 = r(2:2:204);
find(rg1 < thresh)


    L_lochannels = [1:4 55:58 63:66 73 80 82];
    R_lochannels = [51:54 81 89 95:102];
    L_hichannels = [5:25 26 27 28 59:62 67:72 74:75 78];
    R_hichannels = [31:34 36:50 76:77 83:88 90:94];
    mid_channels = [29 30 35 22 79 80];
    all_channels = sort([L_lochannels R_lochannels L_hichannels R_hichannels mid_channels]);
    posUP_chan   = sort([L_lochannels R_hichannels]);
    negUP_chan   = sort([R_lochannels L_hichannels]);

%%
sel1  = ~arrayfun(@(x)(isempty(x.clusters(1).Coh_mean)),SOreference);
sel2  = ~arrayfun(@(x)(isempty(x.clusters(2).Coh_mean)),SOreference);

Cmnegall  = arrayfun(@(x)(x.clusters(1).Coh_mean),SOreference(sel1),'UniformOutput',false);
Cmposall  = arrayfun(@(x)(x.clusters(2).Coh_mean),SOreference(sel2),'UniformOutput',false);
Csnegall  = arrayfun(@(x)(x.clusters(1).Coh_std),SOreference(sel1),'UniformOutput',false);
Csposall  = arrayfun(@(x)(x.clusters(2).Coh_std),SOreference(sel2),'UniformOutput',false);

Cmneg     = cellfun(@mean,Cmnegall);
Cmpos     = cellfun(@mean,Cmposall);
Csneg     = cellfun(@mean,Csnegall);
Cspos     = cellfun(@mean,Csposall);

figure
ch = 1:length(Cmneg);
plot(ch,Cmneg,'b-',ch,Cmpos,'r-'); hline(0,'k')
plot(ch,Cmneg-Cmpos,'*'); hline(0,'k')

figure;
subplot(1,3,1),try hist(Cmneg,40); end
subplot(1,3,2),try hist(Cmpos,40); end
subplot(1,3,3),try hist(Cmneg-Cmpos,40); end


figure;
REF = 1:length(Cmnegall);
nr  = length(REF);
rmax=min(15,nr);
lp  = 1:2:rmax*2;
rp  = 2:2:rmax*2;
for k = 1:nr
  if mod(k,rmax)==1, figure; cnt = 1; end
  ref = REF(k);
  subplot(rmax,2,lp(cnt)),try hist(Cmnegall{ref},40); end; axis tight; title(sprintf('ref %g, neg',ref));
  th = mean(Cmnegall{ref}) + 2*std(Cmnegall{ref}); vline(th,'k');
  Cmnegall{ref} = Cmnegall{ref}(Cmnegall{ref} < th);
  subplot(rmax,2,rp(cnt)),try hist(Cmposall{ref},40); end; axis tight; title('pos')
  th = mean(Cmposall{ref}) + 2*std(Cmposall{ref}); vline(th,'k');
  Cmposall{ref} = Cmposall{ref}(Cmposall{ref} < th);
  cnt = cnt + 1;
end
Cmneg     = cellfun(@mean,Cmnegall(REF));
Cmpos     = cellfun(@mean,Cmposall(REF));
figure
subplot(1,3,1),try hist(Cmneg,40); end
subplot(1,3,2),try hist(Cmpos,40); end
subplot(1,3,3),try hist(Cmneg-Cmpos,40); end


c1 = arrayfun(@(x)(mean(x.clusters(1).Coh_mean)),SOreference);
c2 = arrayfun(@(x)(mean(x.clusters(2).Coh_mean)),SOreference);

%% Plot PSD

% for k = 1:length(SOreference)
%   nn  = length(SOreference(k).clusters(1).epochs);
%   np  = length(SOreference(k).clusters(2).epochs);
%   n   = min(nn,np);
%   nr  = floor(sqrt(n)+1);
%   nc  = nr;
%   Pm1 = zeros(length(SOreference(k).clusters(1).epochs(1).PSDfreq),1); % mean pos
%   Pm2 = zeros(length(SOreference(k).clusters(2).epochs(1).PSDfreq),1); % mean neg
%   figure(1)
%   set(gcf,'Name',SOreference(k).reflabel);
%   for j = 1:n
%     f1  = SOreference(k).clusters(1).epochs(j).PSDfreq;
%     p1  = SOreference(k).clusters(1).epochs(j).PSD;
%     f2  = SOreference(k).clusters(2).epochs(j).PSDfreq;
%     p2  = SOreference(k).clusters(2).epochs(j).PSD;
%     subplot(nr,nc,j)
%     plot(f1,p1,'b-',f2,p2,'r-'); axis tight
%     if j==1,title(sprintf('%g',j)); legend('neg','pos'); end
%     Pm1 = Pm1 + p1/n;
%     Pm2 = Pm2 + p2/n;
%   end
%   figure(2)
%   plot(f1,Pm1,'b-',f2,Pm2,'r-'); axis tight
%   title(sprintf('%s: mean power (n=%g)',SOreference(k).reflabel,n)); legend('neg','pos');
%   pause
% end

% for k = 1:length(SOreference)
%   nr  = 5;
%   nc  = 5;
%   nneg= length(SOreference(k).clusters(1).epochs);
%   npos= length(SOreference(k).clusters(2).epochs);
%   nmax= min(nneg,npos);
%   P0 = zeros(length(SOreference(k).clusters(1).epochs(1).PSDfreq),1);
%   figure(1)
%   set(gcf,'Name',SOreference(k).reflabel);
%   cnt = 0; Pm1=P0; Pm2=P0;
%   for j = 1:nmax
%     cnt = cnt + 1;
%     f1  = SOreference(k).clusters(1).epochs(j).PSDfreq;
%     p1  = SOreference(k).clusters(1).epochs(j).PSD;
%     f2  = SOreference(k).clusters(2).epochs(j).PSDfreq;
%     p2  = SOreference(k).clusters(2).epochs(j).PSD;
%     fix = f1>=foilim(1) & f1<=foilim(2);
%     subplot(nr,nc,cnt)
%     plot(f1,p1,'b-',f2,p2,'r-'); axis tight, title(sprintf('%g: t=%gs',j,t(SOreference(k).peaks(strmatch(SOreference(k).reflabel,{SOreference(k).peaks.label})).negpeak(j))))
%     if cnt==1,legend('neg','pos'); end
%     Pm1 = Pm1 + p1/n;
%     Pm2 = Pm2 + p2/n;
%     if mod(j+1,nr*nc)==1 && j+1<=nmax
%       figure(2)
%       subplot(3,1,1),plot(f1,Pm1,'b-',f2,Pm2,'r-'); axis tight
%       title(sprintf('%s: mean power (n=%g)',SOreference(k).reflabel,n)); legend('neg','pos'); axis tight
%       epochix = j-nr*nc+1:j;
%       nepoch1 = length(epochix);%length(SOreference(k).clusters(1).epochs);
%       Coh1    = SOreference(k).clusters(1).Coh_mean(epochix);
%       PSD1    = SOreference(k).clusters(1).PSD_mean(epochix);
%       nepoch2 = length(epochix);%length(SOreference(k).clusters(2).epochs);
%       Coh2    = SOreference(k).clusters(2).Coh_mean(epochix);
%       PSD2    = SOreference(k).clusters(2).PSD_mean(epochix);
%       subplot(3,1,2),plot(1:nepoch1,Coh1,'b-',1:nepoch2,Coh2,'r-'); title('Coherence'); axis tight
%       subplot(3,1,3),plot(1:nepoch1,PSD1,'b-',1:nepoch2,PSD2,'r-'); title('Power'),xlabel('epoch #'); axis tight
%       pause
%       cnt = 0; Pm1=P0; Pm2=P0;
%       figure(1)
%     end      
%   end
% end
ylim   = [0 1E-19];
foilim = [10 60];
for k = 1:length(SOreference)
  nr  = 5;
  nc  = 5;
  nneg= length(SOreference(k).clusters(1).epochs);
  npos= length(SOreference(k).clusters(2).epochs);
  nmax= min(nneg,npos);
  foi = SOreference(k).clusters(1).epochs(1).PSDfreq;
  fix = foi>=foilim(1) & foi<=foilim(2);
  P0  = zeros(size(SOreference(k).clusters(1).epochs(1).PSDfreq(fix)));
  figure(1)
  set(gcf,'Name',SOreference(k).reflabel);
  cnt = 0; Pm1=P0; Pm2=P0;
  for j = 1:nmax
    cnt = cnt + 1;
    f1  = SOreference(k).clusters(1).epochs(j).PSDfreq(fix);
    p1  = SOreference(k).clusters(1).epochs(j).PSD(fix);
    f2  = SOreference(k).clusters(2).epochs(j).PSDfreq(fix);
    p2  = SOreference(k).clusters(2).epochs(j).PSD(fix);
    subplot(nr,nc,cnt)
    plot(f1,p1,'b-',f2,p2,'r-'); axis tight, set(gca,'ylim',ylim);
    ntk = t(SOreference(k).peaks(strmatch(SOreference(k).reflabel,{SOreference(k).peaks.label})).negpeak(j));
    ptk = t(SOreference(k).peaks(strmatch(SOreference(k).reflabel,{SOreference(k).peaks.label})).pospeak(j));
    title(sprintf('%g: t=%g/%g',j,ntk,ptk))
    if cnt==1,legend('neg','pos'); end
    Pm1 = Pm1 + p1/n;
    Pm2 = Pm2 + p2/n;
    PSDmu1(cnt) = mean(p1);
    PSDmu2(cnt) = mean(p2);
    if mod(j+1,nr*nc)==1 && j+1<=nmax
      figure(2)
      subplot(4,1,1),plot(f1,Pm1,'b-o',f2,Pm2,'r-o'); axis tight, %set(gca,'ylim',ylim);
      title(sprintf('%s: mean power (n=%g)',SOreference(k).reflabel,nr*nc)); legend('neg','pos')
      epochix = j-nr*nc+1:j;
      nepoch1 = length(epochix);%length(SOreference(k).clusters(1).epochs);
%       Coh1    = SOreference(k).clusters(1).Coh_mean(epochix);
      PSD1    = SOreference(k).clusters(1).PSD_mean(epochix);
      nepoch2 = length(epochix);%length(SOreference(k).clusters(2).epochs);
%       Coh2    = SOreference(k).clusters(2).Coh_mean(epochix);
      PSD2    = SOreference(k).clusters(2).PSD_mean(epochix);
%       subplot(4,1,2),plot(1:nepoch1,Coh1,'b-o',1:nepoch2,Coh2,'r-o'); title(sprintf('Coherence (%g-%gHz)',parms.Coh_foilim)); axis tight, set(gca,'ylim',[0 .4]);
      subplot(4,1,3),plot(1:nepoch1,PSD1,'b-o',1:nepoch2,PSD2,'r-o'); title(sprintf('Power (%g-%gHz)',parms.PSD_foilim)),xlabel('epoch #'); axis tight, set(gca,'ylim',ylim);
      subplot(4,1,4),plot(1:nepoch1,PSDmu1,'b-o',1:nepoch2,PSDmu2,'r-o'); title(sprintf('Power (%g-%gHz)',foilim)),xlabel('epoch #'); axis tight, set(gca,'ylim',ylim);
      pause
      cnt = 0; Pm1=P0; Pm2=P0;
      figure(1)
    end      
  end
end




    
    