function [output,info] = p_sa(data,sr,varargin)
% Run analysis on individual time series to extract sympathetic activity
% features (for use on continuous time-series)
% 
% USAGE:
% [output,info] = p_sa(data,sr,TYPE,...)
%
% INPUTS:
% data - NxM matrix of data, with N samples and M subjects
% sr - sample rate of data in Hz
% TYPE - analysis type from following options (and additional inputs):
%
% 'central' - extracts mean (output), as well as median, mode, and standard
% deviation of the (EDA, phasic, tonic, HR) signals
%
% 'auc' - extracts area under the curve of the (EDA, phasic, tonic, HR)
% signals
%
% 'psd',WIN,OVERLAP,F1,F2 - estimates the power spectral density of the
% signals using Welch's periodogram method, using a window of WIN length (in
% seconds) and overlap of OVERLAP seconds. F1 and F2 are the lower and 
% upper frequencies (in Hz) of the frequency range to estimate. If not F1 
% and F2 not specified or both entered as zero, all frequencies are used.
%
% OPTIONAL INPUT: 'progress' - will display progress bar for potentially
% lengthy analysis processes.
%
% OUTPUTS:
% output - measures structure
% info - additional info; analysis parameters etc.

%% PROGRESS BAR CHECK
if ~isempty(find(strcmpi(varargin,'progress'),1)) %use progress bar
    p=1; ix=find(strcmpi(varargin,'progress'),1); varargin(ix)=[];
else
    p=0;
end

%% Centrality measures
if ~isempty(find(strcmpi(varargin,'central'), 1))
    output=nanmean(data);
    info.analysis='central';
    info.std=nanstd(data);
    info.median=nanmedian(data);
    info.mode=mode(data);
    info.samples=length(data(:,1));
    info.missingdata=sum(isnan(data));
end

%% Standard deviation
if ~isempty(find(strcmpi(varargin,'std'), 1))
    output=nanstd(data);
    info.analysis='std';
    info.samples=length(data(:,1));
    info.missingdata=sum(isnan(data));
end

%% Area under Curve
if ~isempty(find(strcmpi(varargin,'auc'),1))
   for ii=1:length(data(1,:)) %for each time-series M
        temp=data(:,ii);
        temp(isnan(temp))=[];
        output(ii)=trapz(temp);
   end 
    info.analysis='auc';
    info.samples=length(data(:,1));
    info.missingdata=sum(isnan(data));
end

%% Power spectral density
if ~isempty(find(strcmpi(varargin,'psd'),1))
    %Checks to ensure correct inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if length(varargin)~=5 && length(varargin)~=3
        error('Incorrect number of inputs for p_sa. Operation cancelled');
    end
    ix=find(strcmpi(varargin,'psd'), 1)+1;
    if length(varargin)==3
        F1=[]; F2=[];
    else
        F1=varargin{ix+2}; F2=varargin{ix+3};
        if F1==0 && F2==0
           F1=[]; F2=[]; 
        end
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
        steps=length(data(1,:));
    end
    info.F=NaN; info.fullF=NaN;
    x=data;
    missing=sum(isnan(x));
    inx=missing>=length(x(:,1))-1;
    x(:,inx)=0;
    x=resample(x,1:length(x(:,1)));
    [pxx,F]=pwelch(x,hann(win),noverlap,win,sr);
    if ~isempty(F1) && ~isempty(F2)
        [~,inx1]=min(abs(F-F1)); %indices of lower freq
        [~,inx2]=min(abs(F-F2)); %indices of higher freq
    else
        inx1=1; inx2=length(F);
    end
    info.totalpower=sum(pxx(inx1:inx2,:),1);
    info.resampled=missing;
    info.samples=length(x(:,1));
    info.F=F(inx1:inx2);
    info.fullF=F;
    info.analysis='psd';
    output=pxx(inx1:inx2,:);
    if p==1, delete(h); end %close waitbar
end