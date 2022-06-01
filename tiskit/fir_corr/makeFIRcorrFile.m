function makeFIRcorrFile(firzeros,sps,filename,timetag)
    % makeFIRcorrFile: create a FIR correction file for Scherbaum routines
    %
    % Usage: makeFIRcorrFile(firzeros,sps,filename,timetag)
    %   firzeros is a list of the FIR zeros
    %   sps is the sampling rate at which the FIR was applied
    %   filename is the name of the desired output file
    %   timetag is the time tag correction in no. of samples
    
    % Find roots of FIR zeros.  Note that these are the INVERSE of the
    % z-plane roots, since the the z-transform uses z^(-n) and "roots" returns
    % the values for z^n.  Therefor the maximum phase coefficients are those
    % INSIDE the unit circle
    fprintf('Calculating roots...'); tic
    z=roots(firzeros);
    
%   UNCOMMENT NEXT TWO LINES TO DEBUG
%    firzeros
%    poly(z)

    fprintf('Done: took %.1f seconds\n',toc);
    lMax= abs(z)<0.999; % Maximum phase is inside unit circle because "roots" gives
                        % roots for positive polynomials (Sherbaum 2007, pp. 126-127)
    lMin= abs(z)>1.001;
    lUC= ~lMax & ~lMin;

    % Plot z-plane zeros
    figure(1);
    zplot(z);
    
    % Plot impulse responses for min Phase, max Phase and unit circle parts
    figure(2);
    impulsePlot(firzeros,z,z(lMax),z(lMin),z(lUC));
    

    % Maximum phase component of filter: NOTE THREE DIFFERENT TRIES BELOW,
    % ONLY LAST TWO WORK
    F_max=z(lMax);       % z transform of maximum filter.  This doesn't duplicate quant40Hz.prt
    F_max=z(lMax).^(-1); % INVERSE z transform of maximum filter: duplicates quant40Hz.prt
    F_max=z(lMin);       % z transform of MINIMUM filter: duplicates quant40Hz.prt
    
    % Calculate equivalent minimum phase filter
    z_mp=z;
    z_mp(lMax)=z(lMax).^(-1);
    fir_mp=fliplr(poly(z_mp));
    fir_mp=fir_mp./(sum(fir_mp))
    figure(3);
         fresp=fft(fir_mp);
         nfft=(length(fresp)+1)/2
       subplot(2,2,1)
            plot(fir_mp,'-o'); axis('tight')
       subplot(2,2,2)
            semilogy(0:1/(nfft-1):1,abs(fresp(1:nfft)),'-')
       subplot(2,2,4)
            plot(0:1/(nfft-1):1,unwrap(angle(fresp(1:nfft))),'-')
  

    fmax=poly(F_max);     % maximum filter, length mx+1
    
    mx=length(F_max);  fprintf('mx:=length(F_max)=%d\n',mx);

    a=-fmax((mx-(1:mx)) +1) ./ fmax((mx) +1); % Eqn 8.13 is missing a neg sign on right side
    b= fmax((0:mx)      +1) ./ fmax((mx) +1);
    
    %% Write the results to a .prt file
    fid=fopen(filename,'w');
    fprintf(fid,'#METHOD: POLYNOMIAL ROOTING\n');
    fprintf(fid,'#FIRNAME file name effective FIR filter coefficients\n');
    fprintf(fid,'test.001\n');
    fprintf(fid,'#NO_EFF: (no. of points for effective FIR filter)\n');
    fprintf(fid,'%d\n',length(firzeros));
    fprintf(fid,'#FDIG_EFF: (digitization frequency in Hz)\n');
    fprintf(fid,'%g\n',sps);
    fprintf(fid,'#FIRDATA_EFF:\n');
    fprintf(fid,'%19.12g %19.12g %19.12g %19.12g %19.12g\n',firzeros);
    fprintf(fid,'\n\n');
    
    fprintf(fid,'#CORRECTION FILTER:\n');
    
    fprintf(fid,'#CORR_AR: (No. of AR coefficients := mx)\n');
    fprintf(fid,'%d\n',mx);
    fprintf(fid,'#CORR_AR_DATA: (AR coefficients a[k] for k = 1 to mx)\n');
    fprintf(fid,'%19.12g %19.12g %19.12g %19.12g %19.12g\n',a);
    fprintf(fid,'\n\n');
    
    fprintf(fid,'#CORR_MA: (No. of MA coefficients := mx+1)\n');
    fprintf(fid,'%d\n',mx+1);
    fprintf(fid,'#CORR_MA_DATA: (MA coefficients b[l] for l = 0 to mx)\n');
    fprintf(fid,'%19.12g %19.12g %19.12g %19.12g %19.12g\n',b);
    fprintf(fid,'\n\n');
    
    fprintf(fid,'#TIMETAG: (time tag correction in no. of samples):\n');
    fprintf(fid,'%g\n',timetag);
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zplot(z)
        
    iMin= find(abs(z)>1.001);
    iMax= find(abs(z)<0.999);
    iCirc=find(abs(z)>=0.999 & abs(z)<=1.001);
    
    ax = gca;
    
    if ~any(imag(z)),
        z = z + 1i*1e-50;
    end;
    newplot;
    if ~isempty(z),
        zh = plot(z,'ko');
        hold on
        plot(z(iMax),'ro');
        plot(z(iMin),'go');
        plot(z(iCirc),'bo');
    else
        zh = [];
    end
    ph = [];
    
    theta = linspace(0,2*pi,70);
    oh = plot(cos(theta),sin(theta),':');
    
    %a = get(ax,'Aspect');
    %set(ax,'Aspect',[a(1),1])
    axis equal
    
    %  zoom out ever so slightly (5%)
    xl=get(ax,'xlim');
    d=diff(xl);
    xl = [xl(1)-.05*d  xl(2)+.05*d];
    set(ax,'xlim',xl);
    yl=get(ax,'ylim');
    d=diff(yl);
    yl = [yl(1)-.05*d  yl(2)+.05*d];
    set(ax,'ylim',yl);
    
    axis(axis);
    
    set(oh,'XData',[get(oh,'XData') NaN ...
        xl(1)-diff(xl)*100 xl(2)+diff(xl)*100 NaN 0 0]);
    set(oh,'YData',[get(oh,'YData') NaN 0 0 NaN ...
        yl(1)-diff(yl)*100 yl(2)+diff(yl)*100]);
    
    handle_counter = 2;
    fuzz = diff(xl)/80; % horiz spacing between 'o' or 'x' and number
    fuzz=0;
    [r,c]=size(z);
    if (r>1)&(c>1),  % multiple columns in z
        ZEE=z;
    else
        ZEE=z(:); c = min(r,c);
    end;
    for which_col = 1:c,      % for each column of ZEE ...
        z = ZEE(:,which_col);
        [mz,z_ind]=mpoles(z);
        for i=2:max(mz),
            j=find(mz==i);
            for k=1:length(j),
                x = real(z(z_ind(j(k)))) + fuzz;
                y = imag(z(z_ind(j(k))));
                if (j(k)~=length(z)),
                    if (mz(j(k)+1)<mz(j(k))),
                        oh(handle_counter) = text(x,y,num2str(i));
                        handle_counter = handle_counter + 1;
                    end
                else
                    oh(handle_counter) = text(x,y,num2str(i));
                    handle_counter = handle_counter + 1;
                end
            end
        end
    end
    set(oh(2:length(oh)),'vertical','bottom');
    
    hold off
    
    xlabel('Real part')
    ylabel('Imaginary part')
end

% Plot impulse responses for min Phase, max Phase and unit circle parts
function impulsePlot(forig,F_total,F_max,F_min,F_UC)
    ftotal=poly(F_total);
    fmax=poly(F_max);
    fmin=poly(F_min);
    fUC=poly(F_UC);

    nx=length(forig);
    subplot(5,1,1); plot(forig,'-o'); title('original FIR response'); set(gca,'XLim',[0 nx]);
    subplot(5,1,2); plot(ftotal,'-o'); title('calculated FIR response');set(gca,'XLim',[0 nx]);
    subplot(5,1,3); plot(fmax,'-o'); title('max phase part');set(gca,'XLim',[0 nx]);
    subplot(5,1,4); plot(fmin,'-o'); title('min phase part');set(gca,'XLim',[0 nx]);
    subplot(5,1,5); plot(fUC,'-o'); title('unit circle part');set(gca,'XLim',[0 nx]);
    
end
