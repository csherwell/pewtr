function [output,info] = p_synchrony(data,sr,varargin)
% Run synchrony analysis (for use on continuous time-series)
% [output,info]= p_synchrony(data,sr,TYPE,____)
%
% INPUTS:
% data - NxM matrix of data, with N samples and M subjects
% sr - sample rate of data in Hz
% TYPE - analysis type from ONE of the following options (and additional
% inputs):
%
% 'SM'- Signal matching. Calculates the absolute pairwise
% differences between each pair of signals and provides the mean distance
% for each pair. Higher values indicate less physiological synchrony.
%
% 'IDM'- Instantaneous derivative matching. Calculates the difference
% between successive samples to find the instanteous derivative (slope) and
% finds the mean of absolute pairwise differences between derivatives for each signal
% pair. Higher values indicate less physiological synchrony.
%
% 'DA' - Directional agreement. Calculates the difference
% between successive samples to find the instanteous derivative (slope) and
% calculates the percentage of matching directions (signs) between each
% pair of signals.
%
% 'Coh',WIN,OVERLAP,F1,F2 - Spectral Coherence Analysis. Uses weighted coherence
% described in Henning et al. (2001). WIN is a numeric input of the Hanning
% window size (in seconds), and OVERLAP is the size (in seconds) of overlap
% for moving window. F1 and F2 are the lower and upper frequencies (in Hz) of
% frequency range to use in coherence estimate. If F1 and F2 not specified,
% all frequencies used.
%
% 'Slope',SWIN,SSTEP,XWIN,XSTEP,MAXLAG(optional),'FZT'(optional) -
% Average Slope Correlation. Calculates the sequential point-by-point differences
% [x(t)-x(t-1)], averages the slope estimates within windows, and cross
% correlates the average slope across subjects. SWIN specifies the size
% (in seconds) of the moving window for calculating average slopes. SSTEP
% specifies the step size (in seconds) of the average slope moving window.
% If SWIN is specified as zero (0), then no windowing is applied and the
% point-by-point differences are calculated between successive samples.
% XWIN specifies the window size for cross correlation, and XSTEP specifies
% the step size for the cross correlation moving window. IF XWIN is
% specified as zero (0), then no moving window is used and the correlation
% is calculated across the entire timeseries.
% OPTIONAL: MAXLAG - Tests the average slope correlation at +/- MAXLAG
% (specified in seconds) to find the maximal correlation and at what time
% lag it occurs.
% OPTIONAL: 'FZT' - Including the argument 'FZT' will perform Fisher's
% Z-transform on correlation estimates.
%
% 'Xcorr',XWIN,XSTEP,MAXLAG(optional),'FZT'(optional) - Cross correlation
% using Pearson correlation coefficient. XWIN specifies the window size for
% moving window cross correlation, and XSTEP specifies the step size for
% the cross correlation moving window. IF XWIN is specified as zero (0),
% then no moving window is used and the correlation is calculated across
% the entire timeseries.
% OPTIONAL: MAXLAG - Tests the cross correlation at +/- MAXLAG
% (specified in seconds) to find the maximal correlation and at what time
% lag it occurs.
% OPTIONAL: 'FZT' - Including the argument 'FZT' will perform Fisher's
% Z-transform on correlation estimates.
%
% 'DTW',MAXLAG - Dynamic time warping. Calculates the minimum sum of
% Euclidean distances between timeseries. MAXLAG specifies, in seconds, the
% constraint that warping is restricted to within MAXLAG. Specifying a
% MAXLAG of zero (0) will apply no constraint (not recommended)
%
% OPTIONAL INPUT: 'progress' - will display progress bar for potentially
% lengthy analysis processes.
%
%
% OUTPUTS:
% output - synchrony measures structure
% info - additional info; analysis parameters etc.

output=[]; info=[];

%% PROGRESS BAR CHECK
if ~isempty(find(strcmpi(varargin,'progress'),1)) %use progress bar
    p=1; ix=find(strcmpi(varargin,'progress'),1); varargin(ix)=[];
else
    p=0;
end

%% SIGNAL MATCHING

if ~isempty(find(strcmpi(varargin,'SM'),1))
    [~,subj]=size(data);
    dist=nan(subj);
    missingdata=nan(subj);
    for ii=1:subj
        for jj=1:subj
            dist(ii,jj)=nanmean(abs(data(:,ii)-data(:,jj)));
            missingdata(ii,jj)=sum(isnan(abs(data(:,ii)-data(:,jj))));
        end
    end
    output=dist;
    info.sr=sr;
    info.samples=length(data(:,1));
    info.missingdata=missingdata;
    info.analysis='SM';
end

%% INSTANTANEOUS DERIVATIVE MATCHING

if ~isempty(find(strcmpi(varargin,'IDM'),1))
    slope=data(2:end,:)-data(1:end-1,:);
    [~,subj]=size(data);
    dist=nan(subj);
    missingdata=nan(subj);
    for ii=1:subj
        for jj=1:subj
            dist(ii,jj)=nanmean(abs(slope(:,ii)-slope(:,jj)));
            missingdata(ii,jj)=sum(isnan(slope(:,ii)-slope(:,jj)));
        end
    end
    output=dist;
    info.sr=sr;
    info.samples=length(slope(:,1));
    info.missingdata=missingdata;
    info.analysis='IDM';
end

%% DIRECTIONAL AGREEMENT

if ~isempty(find(strcmpi(varargin,'DA'),1))
    slope=data(2:end,:)-data(1:end-1,:);
    [~,subj]=size(data);
    da=nan(subj);
    missingdata=nan(subj);
    for ii=1:subj
        for jj=1:subj
            sdir=[slope(:,ii) slope(:,jj)];
            idx=isnan(slope(:,ii)-slope(:,jj));
            sdir(idx,:)=[];
            sdir=sign(sdir);
            da(ii,jj)=sum(sdir(:,1)==sdir(:,2))/length(sdir(:,1));
            missingdata(ii,jj)=sum(isnan(slope(:,ii)-slope(:,jj)));
        end
    end
    output=da;
    info.sr=sr;
    info.samples=length(slope(:,1));
    info.missingdata=missingdata;
    info.analysis='DA';
end

%% SPECTRAL COHERENCE ANALYSIS

if ~isempty(find(strcmpi(varargin,'Coh'), 1))
    %% Checks to ensure correct inputs
    if length(varargin)~=5 && length(varargin)~=3
        error('Incorrect number of inputs for p_synchrony. Operation cancelled');
    end
    
    ix=find(strcmpi(varargin,'Coh'), 1)+1;
    
    if length(varargin)==3
        F1=[]; F2=[];
    else
        F1=varargin{ix+2}; F2=varargin{ix+3};
    end
    
    win=varargin{ix}*sr; %get window size
    if win~=floor(win)
        win=floor(win);
        warning(['Product of window size and samplerate is not an integer. ' ...
            'Adjusted to ' num2str(win/sr) ' seconds']);
    end
    
    noverlap=varargin{ix+1}*sr; %step size
    if noverlap~=floor(noverlap)
        noverlap=floor(noverlap);
        warning(['Product of window size and overlap is not an integer. ' ...
            'Adjusted to ' num2str(win/sr) ' seconds']);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Waitbar
    if p==1
        h=waitbar(0,'1','Name','Running Spectral Coherence');
        subj=length(data(1,:));
        steps=0.5*subj*(subj-1);
        step=0;
    end
    
    %Calculate weighted spectral coherence
    SC=nan(length(data(1,:))); %create output matrix
    info.samples=length(data(:,1));
    info.F=NaN;
    x=data;
    missing=sum(isnan(x));
    inx=missing>=length(x(:,1))-1;
    x(:,inx)=0;
    x=resample(x,1:info.samples);
    [pxx,F]=pwelch(x,hann(win),noverlap,win,sr);
    if ~isempty(F1) && ~isempty(F2)
        [~,inx1]=min(abs(F-F1)); %indices of lower freq
        [~,inx2]=min(abs(F-F2)); %indices of higher freq
    else
        inx1=1; inx2=length(F);
    end
    info.F=F(inx1:inx2);
    pxx=pxx(inx1:inx2,:);
    for ii=1:length(data(1,:)) %for each participant X
        if sum(x(:,ii))==0
            continue
        end
        for kk=1:length(data(1,:)) %for each participant Y
            if isnan(SC(ii,kk)) && ii~=kk
                % Update waitbar and message
                if p==1
                    step=step+1;
                    waitbar(step/steps,h,'Progress')
                end
                [pxy,~]=cpsd(x(:,ii),x(:,kk),hann(win),noverlap,win,sr);
                pxy=pxy(inx1:inx2);
                coh=abs(pxy).^2./(pxx(:,ii).*pxx(:,kk));
                SC(ii,kk)=sum(coh.*(pxx(:,ii).*pxx(:,kk)).^0.5)...
                    /sum((pxx(:,ii).*pxx(:,kk)).^0.5);
                SC(kk,ii)=SC(ii,kk);
            end
        end
    end
    output=SC;
    info.resampled=missing;
    info.analysis='SC';
    if p==1, delete(h); end %close waitbar
end

%% AVERAGE SLOPE ANALYSIS

if ~isempty(find(strcmpi(varargin,'Slope'), 1))
    ix=find(strcmpi(varargin,'Slope'), 1)+1;
    % Adjust slope sample size to integers
    swin=varargin{ix}*sr; %window size in samples
    if swin~=floor(swin)
        swin=floor(swin);
        if swin<1
            warning(['Slope window size too small. Skipping average slope ' ...
                'analysis']);
            return
        else
            warning(['Product of slope window size and samplerate is not an integer. ' ...
                'Adjusted to ' num2str(swin/sr) ' seconds for Average Slope' ...
                ' analysis (Slope sample size)']);
        end
    end
    sstep=varargin{ix+1}*sr; %overlap size in samples
    if sstep~=floor(sstep)
        sstep=floor(sstep);
        if sstep<1
            warning(['Slope step size too small. Skipping average slope ' ...
                'analysis']);
            return
        else
            warning(['Product of slope step size and samplerate is not an integer. ' ...
                'Adjusted to ' num2str(sstep/sr) ' seconds']);
        end
    end
    
    %  Calculate sequential sample-to-sample differences
    slope=NaN(1,length(data(1,:))); %first row empty
    slope=[slope; data(2:end,:)-data(1:end-1,:)];
    if swin~=0
        steps=1:sstep:length(slope(:,1))-swin+1;
        temp=NaN(length(steps),length(slope(1,:)));
        for kk=1:length(steps)
            temp(kk,:)=nanmean(slope(steps(kk):steps(kk)+sstep-1,:));
        end
        slope=temp;
        sr=sr/sstep; %change sample rate to reflect stepped sampling
    end
    
    % Adjust correlation sample size
    xwin=varargin{ix+2}*sr; %window size in samples
    if xwin~=floor(xwin)
        xwin=floor(xwin);
        if xwin<1
            error(['Correlation window size too small (less than one ' ...
                'average slope sample). Operation cancelled']);
        else
            warning(['Product of correlation window size and samplerate is ' ...
                'not an integer. Adjusted to ' num2str(xwin/sr) ' seconds']);
        end
    end
    xstep=varargin{ix+3}*sr; %step size in samples
    if xstep~=floor(xstep)
        xstep=floor(xstep);
        if xstep<1
            error(['Correlation step size too small (less than one ' ...
                'average slope sample). Operation cancelled']);
        else
            warning(['Product of correlation step size and samplerate is ' ...
                'not an integer. Adjusted to ' num2str(xstep/sr) ' seconds']);
        end
    end
    maxlag=0;
    if length(varargin)>5
        if isnumeric(varargin{ix+4})
            maxlag=varargin{ix+4}*sr;
            if maxlag~=floor(maxlag)
                maxlag=floor(maxlag);
                warning(['Product of lag and samplerate is not an integer. ' ...
                    'Adjusted to ' num2str(maxlag/sr) ' seconds for Average' ...
                    ' Slope Analysis (lagged correlations)']);
            end
        end
    end
    if ~isempty(find(strcmpi(varargin,'FZT'), 1)) %check for fishers transform
        fzt=1;
    else
        fzt=0;
    end
    
    %Cross correlate (Zero-lag)
    subj=length(slope(1,:));
    if maxlag==0
        if xwin==0 %full window correlation
             h=waitbar(0,'1','Name',['Running Average Slope Correlation' ...
                    '(zero-lag)']);
            rval(:,:,1)=corr(slope,'rows','pairwise');
            xstep=1;
        else
            samples=length(1:xstep:length(slope(:,1))-xwin+1);
            rval=NaN(subj,subj,samples);
            ss=1; %steps
            
            if p==1 % Waitbar
                h=waitbar(0,'1','Name',['Running Average Slope Correlation' ...
                    '(zero-lag)']);
            end
            
            for kk=1:samples
                if p==1, waitbar(kk/samples,h,'Progress'); end %waitbar update
                rval(:,:,kk)=corr(slope(ss:ss+xwin-1,:),'rows','pairwise');
                ss=ss+xstep;
            end
        end
        if fzt==1, rval=atanh(rval); end %Fishers z-transform
        
    elseif maxlag~=0 %Lagged correlations
        slope=cat(1,zeros(maxlag,subj),slope);
        slope=cat(1,slope,zeros(maxlag,subj));
        if xwin==0 %full window correlation
            rval=NaN(subj,subj,1);
            delay=rval;
            xstep=1;
            if p==1 %waitbar
                h=waitbar(0,'1','Name','Running Cross Correlation');
                steps=((subj*subj)-subj)*((maxlag*2)+1);
                step=0;
            end
            for ii=1:subj
                for jj=1:subj
                    temp=nan(1,(maxlag*2)+1);
                    for lag=-maxlag:maxlag
                        if p==1, step=step+1; waitbar(step/steps,h,'Progress'); end %waitbar
                        temp(lag+maxlag+1)=corr(slope(lag+maxlag+1:end+lag-maxlag,ii),...
                            slope(maxlag+1:end-maxlag,jj),'rows','pairwise');
                    end
                    if fzt==1, temp=atanh(temp); end %Fishers z-transform
                    rval(ii,jj,1)=max(temp); %log maximum R val
                    idx=find(temp==max(temp)); %index of max R val
                    if ~isempty(idx)
                        [~,idi]=min(abs(idx-maxlag-1)); %if >1, find smallest delay
                        delay(ii,jj,1)=idx(idi)-maxlag-1;
                    end
                end
            end
        else %Windowed correlation
            samples=maxlag+1:xstep:length(slope(:,1))-xwin-maxlag;
            rval=NaN(subj,subj,length(samples));
            delay=NaN(subj,subj,length(samples));
            if p==1 %waitbar
                h=waitbar(0,'1','Name',['Running Average Slope Correlation' ...
                    '(' num2str(maxlag/sr) 'sec max lag)']);
                steps=((subj*subj)-subj)*((maxlag*2)+1)*length(samples);
                step=0;
            end
            for ii=1:subj
                for jj=1:subj
                    if ii~=jj
                        ss=1; %steps
                        for kk=samples %for each xstep
                            temp=nan(1,length(-maxlag:maxlag));
                            for lag=-maxlag:maxlag
                                if p==1, step=step+1; waitbar(step/steps,h,'Progress'); end %waitbar
                                temp(lag+maxlag+1)=corr(slope(kk+lag:kk+lag+xwin-1,ii),...
                                    slope(kk:kk+xwin-1,jj),'rows','pairwise');
                            end
                            if fzt==1, temp=atanh(temp); end %Fishers transform
                            rval(ii,jj,ss)=max(temp);
                            if ~isnan(max(temp))
                            idx=find(temp==max(temp));
                            [~,idi]=min(abs(idx-maxlag-1));
                            delay(ii,jj,ss)=idx(idi)-maxlag-1;
                            end
                            ss=ss+1;
                        end
                    end
                end
            end
        end
    end
    
    %pairwise missing data
    missingdata=nan(subj);
    for ii=1:subj
        for jj=1:subj
            missingdata(ii,jj)=sum(isnan(data(:,ii)-data(:,jj)));
        end
    end
    
    output=rval;
    info.sr=sr/xstep; %change sample rate to reflect stepped sampling
    info.maxlag=maxlag;
    if maxlag~=0, info.delay=delay; end
    info.samples=length(data(:,1));
    info.missingdata=missingdata;
    info.xwin=varargin{ix+2};
    info.xstep=varargin{ix+3};
    info.swin=varargin{ix};
    info.sstep=varargin{ix+1};
    info.analysis='slope';
    info.fzt=fzt;
    
    if p==1, delete(h); end %close waitbar
end


%% CROSS CORRELATION ANALYSIS

if ~isempty(find(strcmpi(varargin,'Xcorr'), 1))
    ix=find(strcmpi(varargin,'Xcorr'), 1)+1;
    % Adjust windows to integers
    xwin=varargin{ix}*sr; %window size in samples
    if xwin~=floor(xwin)
        xwin=floor(xwin);
        if xwin<1
            error(['Correlation window size too small (less than one ' ...
                'sample). Operation cancelled']);
        else
            warning(['Product of window size and samplerate is not an integer. ' ...
                'Adjusted to ' num2str(xwin/sr) ' seconds']);
        end
    end
    xstep=varargin{ix+1}*sr; %overlap size in samples
    if xstep~=floor(xstep)
        xstep=floor(xstep);
        if xstep<1
            error(['Correlation step size too small (less than one ' ...
                'sample). Operation cancelled']);
        else
            warning(['Product of step size and samplerate is not an integer. ' ...
                'Adjusted to ' num2str(xstep/sr) ' seconds']);
        end
    end
    maxlag=0;
    if length(varargin)>3 %lagged correlations
        if isnumeric(varargin{ix+2})
            maxlag=varargin{ix+2}*sr;
            if maxlag~=floor(maxlag)
                maxlag=floor(maxlag);
                warning(['Product of lag and samplerate is not an integer. ' ...
                    'Adjusted to ' num2str(maxlag/sr) ' seconds']);
            end
        end
    end
    if ~isempty(find(strcmpi(varargin,'FZT'), 1)) %check for fishers transform
        fzt=1;
    else
        fzt=0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Cross correlate (Zero-lag)
    subj=length(data(1,:));
    if maxlag==0
        if xwin==0 %full window correlation
            if p==1 % Waitbar
                h=waitbar(0,'1','Name',['Running Cross Correlation' ...
                    '(zero-lag)']);
            end
            rval(:,:,1)=corr(data,'rows','pairwise');
            xstep=1;
        else
            samples=length(1:xstep:length(data(:,1))-xwin+1);
            rval=NaN(subj,subj,samples);
            ss=1; %steps
            
            if p==1 % Waitbar
                h=waitbar(0,'1','Name',['Running Cross Correlation' ...
                    '(zero-lag)']);
            end
            
            for kk=1:samples
                if p==1, waitbar(kk/samples,h,'Progress'); end %waitbar update
                rval(:,:,kk)=corr(data(ss:ss+xwin,:),'rows','pairwise');
                ss=ss+xstep;
            end
        end
        if fzt==1, rval=atanh(rval); end %Fishers z-transform
        %Lagged correlations
    elseif maxlag~=0
        data=cat(1,zeros(maxlag,subj),data);
        data=cat(1,data,zeros(maxlag,subj));
        if xwin==0 %full window correlation
            rval=NaN(subj,subj,1);
            delay=rval;
            xstep=1;
            if p==1 %waitbar
                h=waitbar(0,'1','Name','Running Cross Correlation');
                steps=((subj*subj)-subj)*((maxlag*2)+1);
                step=0;
            end
            for ii=1:subj
                for jj=1:subj
                    temp=nan(1,(maxlag*2)+1);
                    for lag=-maxlag:maxlag
                        if p==1, step=step+1; waitbar(step/steps,h,'Progress'); end %waitbar
                        temp(lag+maxlag+1)=corr(data(lag+maxlag+1:end+lag-maxlag,ii),...
                            data(maxlag+1:end-maxlag,jj),'rows','pairwise');
                    end
                    if fzt==1, temp=atanh(temp); end %Fishers z-transform
                    rval(ii,jj,1)=max(temp); %log maximum R val
                    idx=find(temp==max(temp)); %index of max R val
                    if ~isempty(idx)
                        [~,idi]=min(abs(idx-maxlag-1)); %if >1, find smallest delay
                        delay(ii,jj,1)=idx(idi)-maxlag-1;
                    end
                end
            end
        else %Windowed correlation
            samples=maxlag+1:xstep:length(data(:,1))-xwin-maxlag;
            rval=NaN(subj,subj,length(samples));
            delay=rval;
            if p==1 %waitbar
                h=waitbar(0,'1','Name',['Running Cross Correlation' ...
                    '(' num2str(maxlag/sr) 'sec max lag)']);
                steps=((subj*subj)-subj)*((maxlag*2)+1)*length(samples);
                step=0;
            end
            for ii=1:subj
                for jj=1:subj
                    ss=1; %steps
                    for kk=samples %for each xstep
                        temp=nan(1,(maxlag*2)+1);
                        for lag=-maxlag:maxlag
                            if p==1, step=step+1; waitbar(step/steps,h,'Progress'); end %waitbar
                            temp(lag+maxlag+1)=corr(data(kk+lag:kk+lag+xwin-1,ii),...
                                data(kk:kk+xwin-1,jj),'rows','pairwise');
                        end
                        if fzt==1, temp=atanh(temp); end %Fishers z-transform
                        rval(ii,jj,ss)=max(temp); %log maximum R val
                        idx=find(temp==max(temp)); %index of max R val
                        if ~isempty(idx)
                            [~,idi]=min(abs(idx-maxlag-1)); %if >1, find smallest delay
                            delay(ii,jj,ss)=idx(idi)-maxlag-1;
                        end
                        ss=ss+1;
                    end
                end
            end
        end
    end
    
    %pairwise missing data
    missingdata=nan(subj);
    for ii=1:subj
        for jj=1:subj
            missingdata(ii,jj)=sum(isnan(data(:,ii)-data(:,jj)));
        end
    end
    
    output=rval;
    info.sr=sr/xstep; %change sample rate to reflect stepped sampling
    info.maxlag=maxlag;
    if maxlag~=0, info.delay=delay; end
    info.samples=length(data(:,1));
    info.missingdata=missingdata;
    info.xwin=varargin{ix};
    info.xstep=varargin{ix+1};
    info.analysis='xcorr';
    info.fzt=fzt;
    
    if p==1, delete(h); end %close waitbar
end

%% DYNAMIC TIME WARPING

if ~isempty(find(strcmpi(varargin,'dtw'), 1))
    % Checks to ensure correct inputs
    if length(varargin)~=2
        error('Incorrect number of inputs for p_synchrony. Operation cancelled');
    end
    ix=find(strcmpi(varargin,'dtw'), 1)+1;
    if ~isnumeric(varargin{ix})
        error('Incorrect inputs for p_synchrony. Operation cancelled');
    end
    maxlag=varargin{ix}*sr; %window size in samples
    if maxlag~=floor(maxlag)
        maxlag=floor(maxlag);
        warning(['Product of lag constraint and samplerate is not an integer. ' ...
            'Adjusted to ' num2str(maxlag/sr) ' seconds']);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %run DTW
    subj=length(data(1,:));
    dist=nan(subj);
    % Waitbar
    if p==1
        h=waitbar(0,'1','Name','Running Dynamic Time Warping');
        steps=subj*((subj+1)/2)-subj;
        step=0;
    end
    for ii=1:subj
        X=data(:,ii); X(isnan(X))=0;
        for jj=1:subj
            if ii~=jj && isnan(dist(ii,jj))
                if p==1 %waitbar
                    step=step+1;
                    waitbar(step/steps,h,'Progress')
                end
                Y=data(:,jj); Y(isnan(Y))=0;
                if maxlag==0
                    dist(ii,jj)=dtw(X,Y);
                    dist(jj,ii)=dist(ii,jj);
                else
                    dist(ii,jj)=dtw(X,Y,maxlag);
                    dist(jj,ii)=dist(ii,jj);
                end
            end
        end
    end
    
    output=dist;
    info.maxlag=maxlag;
    info.samples=length(data(:,1));
    info.missingdata=sum(isnan(data));
    info.analysis='dtw';
    if p==1, delete(h); end %close waitbar
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(output)
    error('Incorrect input arguments. No analyses performed')
end

