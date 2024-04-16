function [pewtr4] = p_analysis(pewtr3,datatype,conditions,varargin)
% Batch function for running synchrony / sympathetic activity analyses on
% sochro format data structures.
% [sochro4] = sochro_analysis(sochro, datatype,conditions,________)
%
% INPUTS:
% sochro - sochro3 (preprocessed) format data structure
%
% datatype - 'EDA' (electrodermal activity), 'PHA' (phasic activity),
% 'TON' (tonic activity), or 'HR' (heartrate) to be analysed. For multiple
% datatypes, enter as cell array: {'EDA', 'PHA', 'TON'}
%
% conditions - numerical array of conditions indices to be analyzed -
% leaving empty will run analyses on each condition (including the entire
% time-series). Include zero [0] in an array of indices to include the
% entire time-series as an N+1 condition for N conditions. Leave zero [0]
% out of an array of indices to skip whole time-series analysis.
%
% OPTIONS:
%
% - SYNCHRONY ANALYSES -
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
% frequency range to use in coherence estimate. Specify F1 and F2 as
% zero(0) to use all frequencies.
%
% 'Slope',SWIN,SSTEP,XWIN,XSTEP,MAXLAG(optional) - Average Slope
% Correlation. Calculates the sequential point-by-point differences
% [x(t)-x(t-1)], averages the slope estimates within windows, and cross
% correlates the average slope across subjects. SWIN specifies the size
% (in seconds) of the moving window for calculating average slopes. SSTEP
% specifies the step size (in seconds) of the average slope moving window.
% If SWIN is specified as zero (0), then no windowing is applied and the
% point-by-point differences are calculated between successive samples.
% XWIN specifies the window size for moving window cross correlation,
% and XSTEP specifies the step size for the cross correlation moving
% window. IF XWIN is specified as zero (0), then no moving window is used
% and the correlation is calculated across the entire timeseries.
% OPTIONAL: MAXLAG - Tests the average slope correlation at +/- MAXLAG
% (specified in seconds) to find the maximal correlation and at what time
% lag it occurs.
% OPTIONAL - replace 'Slope' with 'FSlope' to perform Fisher's z-transform
% on the resulting correlation coefficients.
%
% 'Xcorr',XWIN,XSTEP,MAXLAG(optional) - Cross correlation using Pearson
% correlation coefficients. XWIN specifies time windows (in seconds) with
% which to calculate cross-correlations.
% Default is 0, where no windowing is applied and cross correlation is
% calculated over entire time series. XSTEP determines the step size
% (in seconds) of the moving correlation window.
% OPTIONAL: MAXLAG - Tests the cross correlation at +/- MAXLAG
% (specified in seconds) to find the maximal correlation and at what time
% lag it occurs.
% OPTIONAL - replace 'Xcorr' with 'FXcorr' to perform Fisher's z-transform
% on the resulting correlation coefficients.
%
% 'DTW',MAXLAG - Dynamic time warping. Calculates the minimum sum of
% Euclidean distances between timeseries. MAXLAG specifies, in seconds, the
% constraint that warping is restricted to within MAXLAG. Specifying a
% MAXLAG of zero (0) will apply no constraint (not recommended)
%
% - SYMPATHETIC ACTIVITY ANALYSES -
%
% 'Central' - extracts mean, median, mode, and standard deviation of the 
% (EDA, phasic, tonic, HR) signals
%
% 'AUC' - extracts area under the curve of the (EDA, phasic, tonic, HR)
% signals
%
% 'PSD',WIN,OVERLAP,F1,F2 - estimates the power spectral density of the
% signals using Welch's periodogram method, using a window of WIN length (in
% seconds) and overlap of OVERLAP seconds. F1 and F2 are the lower and
% upper frequencies (in Hz) of the frequency range to estimate. If inputs
% for F1 and F2 are zero (0), the entire spectrum is used.
%
% 'SCR',WIN,STEP - extracts skin conductance response measures from phasic component
% of the EDA signal. Phasic data must be present to extracts SCR measures.
% WIN specifies a moving window (of WIN length in seconds) with which to
% calculate mean peaks per window, with a step size of STEP (in
% seconds). Using a WIN size of zero (0) will calculate mean as number of
% peaks divided by data length.
%
% - OTHER OPTIONS -
%
% 'progress' - displays pop-up progress bar

%% INPUT CHECK
if isempty(pewtr3)
    error('Missing input for p_analysis: missing pewtr3 data structure')
elseif isempty(datatype)
    error('Missing input for p_analysis: no datatypes specified')
elseif isempty(find(strcmpi(datatype,'EDA'), 1)) && ...
        isempty(find(strcmpi(datatype,'PHA'), 1)) && ...
        isempty(find(strcmpi(datatype,'TON'), 1)) && ...
        isempty(find(strcmpi(datatype,'HR'), 1))
    error('Incorrect input for p_analysis: incorrect datatype specified')
elseif isempty(varargin)
    error('Incorrect input for p_analysis: no analyses specified')
end

if isempty(conditions)
    cinx=[0 2:length(pewtr3.conditions)];
else
    cinx=conditions;
end

if ~iscell(datatype)
    datatype={datatype};
end

if ~isempty(find(strcmpi(varargin,'progress'),1)) %use progress bar
    p=1; ix=find(strcmpi(varargin,'progress'),1); varargin(ix)=[];
else
    p=0;
end
analyses=find(cellfun(@(c) ischar(c),varargin)); %inx of analysis args
output = struct('datatype',datatype,'analysis',[]);

% Waitbar
if p==1
    fname=strsplit(pewtr3.file.sdir,{'/' '\'}); fname=fname{end};
    f=waitbar(0,'1','Name',['Analyzing dataset: ' fname ]);
    steps=length(analyses)*length(datatype)*length(cinx);
    step=0;
end


for dd=1:length(datatype) %for each datatype
    disp(['PEWTR Analysis: running '  datatype{dd} ' data']);
    
    % Update waitbar and message
    if p==1
        waitbar(step/steps,f,['running '  datatype{dd} ' data'])
    end
    
    sigpro=[]; %for sigprocessing history
    %% ASSIGN DATA
    switch datatype{dd}
        case {'EDA','eda'}
            data=pewtr3.data.EDA;
            sr=pewtr3.file.Hz.EDAhz;
            time=pewtr3.time.EDAmins;
        case {'HR','hr'}
            data=pewtr3.data.HR;
            sr=pewtr3.file.Hz.HRhz;
            time=pewtr3.time.HRmins;
        case {'PHA','pha'}
            if ~isfield(pewtr3.data,'Phasic')
                error('Error in p_analysis: no phasic data present');
            end
            data=pewtr3.data.Phasic;
            sr=pewtr3.file.Hz.Phasichz;
            time=pewtr3.time.Phasicmins;
        case {'TON','ton'}
            if ~isfield(pewtr3.data,'Tonic')
                error('Error in p_analysis: no tonic data present');
            end
            data=pewtr3.data.Tonic;
            sr=pewtr3.file.Hz.Tonichz;
            time=pewtr3.time.Tonicmins;
    end
    subj=length(data(1,:));
    
    %% GET CONDITION WINDOWS
    %Condition window indices and label (tinx rows: indices; label;
    % minutes). Save to output structure
    tinx=cell(2,length(cinx));
    for kk=1:length(cinx) %for each condition
        if cinx(kk)==0
            tinx{1,kk}={[1, length(time)]};
            tinx{2,kk}='Session';
            tinx{3,kk}={[time(1) time(end)]};
            output(dd).condition(kk).label='Session';
            output(dd).condition(kk).inx={[1, length(time)]};
            output(dd).condition(kk).t={[time(1) time(end)]};
            output(dd).condition(kk).data{1}=data;
        else
            tinx(2,kk)=pewtr3.conditions(cinx(kk),1);
            output(dd).condition(kk).label=tinx{2,kk};
            for ww=2:length(pewtr3.conditions(cinx(kk),:)) %for each window
                [~,tinx{1,kk}{ww-1}(1)]=min(abs(minutes(time)- ...
                    pewtr3.conditions{cinx(kk),ww}(1)));
                [~,tinx{1,kk}{ww-1}(2)]=min(abs(minutes(time)- ...
                    pewtr3.conditions{cinx(kk),ww}(2)));
                tinx{3,kk}{ww}=pewtr3.conditions(cinx(kk),ww);
                output(dd).condition(kk).inx{ww}=tinx{1,kk}{ww-1};
                output(dd).condition(kk).t{ww}=tinx{3,kk}{ww};
                output(dd).condition(kk).data{ww}=data(tinx{1,kk}{ww-1}(1):tinx{1,kk}{ww-1}(2),:);
            end
        end
    end
    
    %% GET ANALYSIS ARGUMENTS
    for aa = 1:length(analyses)  %for each analysis type
        output(dd).analysis(aa).parameters=[];
        switch lower(varargin{analyses(aa)})
            case 'sm'
                amsg='Running Signal Matching analysis';
                disp(amsg)
                output(dd).analysis(aa).label='sm';
                args={'sm'}; type='s';
            case 'idm'
                amsg=('Running Instantaneous Derivative Matching analysis');
                disp(amsg)
                output(dd).analysis(aa).label='idm';
                args={'idm'}; type='s';
            case 'da'
                amsg=('Running Directional Agreement analysis');
                disp(amsg)
                output(dd).analysis(aa).label='da';
                args={'da'}; type='s';
            case 'coh'
                amsg=('Running Spectral Coherence analysis');
                disp(amsg)
                output(dd).analysis(aa).label='coh';
                ix=find(strcmpi(varargin,'coh'), 1); %get analysis arguments
                output(dd).analysis(aa).parameters= struct('win', ...
                    varargin{ix+1},'overlap',varargin{ix+2},'F1', ...
                    varargin{ix+3},'F2',varargin{ix+4});
                args=varargin(ix:ix+4); type='s';
            case {'slope','fslope'}
                amsg=('Running Average Slope analysis');
                disp(amsg)
                ix=find(strcmpi(varargin,'slope'), 1); %get analysis arguments
                output(dd).analysis(aa).label='slope';
                fzt=0;
                args={'slope'};
                if isempty(ix)
                    ix=find(strcmpi(varargin,'fslope'), 1);
                    fzt=1;
                end
                if length(varargin)<ix+5
                    output(dd).analysis(aa).parameters= struct('swin', ...
                        varargin{ix+1},'sstep',varargin{ix+2},'xwin', ...
                        varargin{ix+3},'xstep',varargin{ix+4},'maxlag',0,...
                        'fzt',fzt);
                    args=[args varargin(ix+1:ix+4)];
                elseif ~isnumeric(varargin{ix+5})
                    output(dd).analysis(aa).parameters= struct('swin', ...
                        varargin{ix+1},'sstep',varargin{ix+2},'xwin', ...
                        varargin{ix+3},'xstep',varargin{ix+4},'maxlag',0,...
                        'fzt',fzt);
                    args=[args varargin(ix+1:ix+4)];
                else
                    output(dd).analysis(aa).parameters= struct('swin', ...
                        varargin{ix+1},'sstep',varargin{ix+2},'xwin', ...
                        varargin{ix+3},'xstep',varargin{ix+4},'maxlag',...
                        varargin{ix+5},'fzt',fzt);
                    args=[args varargin(ix+1:ix+5)];
                end
                if fzt==1, args{end+1}='fzt'; end
                type='s';
            case {'xcorr','fxcorr'}
                amsg=('Running Cross Correlation analysis');
                disp(amsg)
                ix=find(strcmpi(varargin,'xcorr'), 1); %get analysis arguments
                output(dd).analysis(aa).label='xcorr';
                fzt=0;
                args={'xcorr'};
                if isempty(ix)
                    ix=find(strcmpi(varargin,'fxcorr'), 1);
                    fzt=1;
                end
                if length(varargin)<ix+3
                    output(dd).analysis(aa).parameters= struct('xwin', ...
                        varargin{ix+1},'xstep',varargin{ix+2},'maxlag',0,...
                        'fzt',fzt);
                    args=[args varargin(ix+1:ix+2)];
                elseif ~isnumeric(varargin{ix+3})
                    output(dd).analysis(aa).parameters= struct('xwin', ...
                        varargin{ix+1},'xstep',varargin{ix+2},'maxlag',0,...
                        'fzt',fzt);
                    args=[args varargin(ix+1:ix+2)];
                else
                    output(dd).analysis(aa).parameters= struct('xwin', ...
                        varargin{ix+1},'xstep',varargin{ix+2},'maxlag',...
                        varargin{ix+3},'fzt',fzt);
                    args=[args varargin(ix+1:ix+3)];
                end
                if fzt==1, args{end+1}='fzt'; end
                type='s';
            case 'dtw'
                amsg=('Running Dynamic time warping analysis');
                disp(amsg)
                output(dd).analysis(aa).label='dtw';
                ix=find(strcmpi(varargin,'dtw'), 1); %get analysis arguments
                output(dd).analysis(aa).parameters= struct('maxlag', ...
                    varargin{ix+1});
                args=varargin(ix:ix+1); type='s';
            case 'central'
                amsg=('Extracting centrality measures of signal');
                disp(amsg)
                output(dd).analysis(aa).label='central';
                args={'central'}; type='i';
            case 'auc'
                amsg=('Running Area under curve estimation');
                disp(amsg)
                output(dd).analysis(aa).label='auc';
                args={'auc'}; type='i';
            case 'psd'
                amsg=('Running Power spectral density estimation');
                disp(amsg)
                output(dd).analysis(aa).label='psd';
                ix=find(strcmpi(varargin,'psd'), 1);
                output(dd).analysis(aa).parameters= struct('win', ...
                    varargin{ix+1},'overlap',varargin{ix+2},'F1', ...
                    varargin{ix+3},'F2',varargin{ix+4});
                args=varargin(ix:ix+4); type='i';
            case 'scr'
                if strcmpi(datatype{dd},'PHA')
                    amsg=('Extracting SCR measures');
                    disp(amsg)
                    output(dd).analysis(aa).label='scr';
                    ix=find(strcmpi(varargin,'scr'), 1);
                    output(dd).analysis(aa).parameters= struct('win', ...
                        varargin{ix+1},'step',varargin{ix+2});
                    args=varargin(ix+1:ix+2); type='scr';
                else
                    continue
                end
        end
        output(dd).analysis(aa).timestamp=datestr(now);
        output(dd).analysis(aa).condition=struct;
        if p==1, args{end+1}='progress'; end
        
        %% RUN ANALYSIS PER WINDOW, PER CONDITION
        for kk=1:length(cinx) %for each condition
%             try
                temp=struct;
                for ww=1:length(tinx{1,kk}) %for each window
                    % Update waitbar and message
                    if p==1
                        step=step+(1/length(tinx{1,kk}));
                        waitbar(step/steps,f,[amsg ': condition ' num2str(kk) ...
                            ' of ' num2str(length(cinx)) ]);
                    end
                    x=data(tinx{1,kk}{ww}(1):tinx{1,kk}{ww}(2),:); %get window data
                    if strcmpi(type,'s') %Synchrony analysis
                        [temp(ww).est,temp(ww).tinfo]=p_synchrony(x,sr,args{:});
                    elseif strcmpi(type,'i') %Individual sympathetic activity
                        [temp(ww).est,temp(ww).tinfo]=p_sa(x,sr,args{:});
                    elseif strcmpi(type,'scr') %Skin conductance responses
                        [temp(ww).est,temp(ww).tinfo]=p_scr(pewtr3.SCR,...
                            tinx{1,kk}{ww},args{1:2});
                        temp(ww).est.samples=length(x(:,1));
                        temp(ww).est.missingdata=sum(isnan(x));
                    end
                end %end of window loop
                
                %% Combine windows and add data to output structure
                windows=struct; tinfo=[temp.tinfo];
                for ww=1:length(tinx{1,kk})
                    windows(ww).inx=tinx{1,kk}{ww};
                    windows(ww).t=tinx{3,kk}{ww};
                end
                
                switch lower(varargin{analyses(aa)})
                    case {'sm','idm','da'}
                        w=cat(3,tinfo.samples)-cat(3,tinfo.missingdata); %samples
                        C=cat(3,temp.est); %sm per window
                        for ww=1:length(tinx{1,kk})
                            windows(ww).time=time(tinx{1,kk}{ww}(1):tinx{1,kk}{ww}(2));
                            windows(ww).est=temp(ww).est;
                            windows(ww).samples=tinfo(ww).samples;
                            windows(ww).missingdata=tinfo(ww).missingdata;
                        end
                        stemp=struct('label',tinx(2,kk),...
                            'est',sum(C.*w,3)./sum(w,3), ... %weighted mean
                            'samples', sum([tinfo.samples]), ...
                            'missingdata', sum(cat(3,tinfo.missingdata),3), ...
                            'windows', windows);
                        
                    case 'coh'
                        w=[];
                        w(1,1,:)=[tinfo.samples]; %samples
                        C=cat(3,temp.est); %coherence per window
                        for ww=1:length(tinx{1,kk})
                            windows(ww).time=time(tinx{1,kk}{ww}(1):tinx{1,kk}{ww}(2));
                            windows(ww).est=temp(ww).est;
                            windows(ww).samples=tinfo(ww).samples;
                            windows(ww).resampled=tinfo(ww).resampled;
                        end
                        stemp=struct('label',tinx(2,kk),...
                            'est',sum(C.*w,3)./sum(w,3), ... %weighted mean
                            'samples', sum(w,3), ...
                            'resampled', sum(cat(1,tinfo.resampled),1), ...
                            'F', tinfo(1).F, ...
                            'windows', windows);
                        
                    case {'xcorr','fxcorr','slope','fslope'}
                        if sum([tinfo.xwin])==0 %no windowing applied
                            w=cat(3,tinfo.samples)-cat(3,tinfo.missingdata); %samples
                            C=cat(3,temp.est);
                            for ww=1:length(tinx{1,kk})
                                windows(ww).est=temp(ww).est;
                                windows(ww).samples=tinfo(ww).samples;
                                windows(ww).missingdata=tinfo(ww).missingdata;
                                if tinfo(1).maxlag~=0
                                   windows(ww).delay=tinfo(ww).delay; 
                                end
                            end
                            stemp=struct('label',tinx(2,kk),...
                                'tseries',[], ...
                                'timeline',[],...
                                'mean',sum(C.*w,3)./sum(w,3), ... %weighted mean
                                'std',[],...
                                'sr', tinfo(1).sr, ...
                                'samples', sum([tinfo.samples]), ...
                                'missingdata', sum(cat(3,tinfo.missingdata),3), ...
                                'windows', windows);
                            if tinfo(1).maxlag~=0
                                C=cat(3,tinfo.delay);
                                stemp.delay=sum(C.*w,3)./sum(w,3);
                            end
                        else %windowed correlation
                            timeline=downsample(time,sr/tinfo(1).sr);
                            tseries=nan(subj,subj,length(timeline));
                            for ww=1:length(tinx{1,kk}) %for each window
                                [~,t1]=min(abs(timeline-tinx{1,kk}{ww}(1))); %start idx
                                tseries(:,:,t1:t1+length(temp(ww).est(1,1,:))-1) = ...
                                    temp(ww).est;
                                windows(ww).est=temp(ww).est;
                                windows(ww).samples=tinfo(ww).samples;
                                windows(ww).missingdata=tinfo(ww).missingdata;
                                if tinfo(1).maxlag~=0
                                   windows(ww).delay=tinfo(ww).delay; 
                                end
                            end
                            %SSI calculation
                            R1=cat(3,temp.est); R2=R1;
                            R1(R1<0)=0; R1=nansum(R1,3);
                            R2(R2>0)=0; R2=abs(nansum(R2,3));
                            SSI=log(R1./R2);
                            stemp=struct('label',tinx(2,kk),...
                                'tseries',tseries, ...
                                'timeline',timeline,...
                                'mean',nanmean(cat(3,temp.est),3), ...
                                'median',nanmedian(cat(3,temp.est),3),...
                                'std',nanstd(cat(3,temp.est),0,3),...
                                'SSI',SSI,...
                                'sr', tinfo(1).sr, ...
                                'samples', sum([tinfo.samples]), ...
                                'missingdata', sum(isnan(cat(3,temp.est)),3), ...
                                'windows', windows);
                            if tinfo(1).maxlag~=0
                                stemp.delay=nanmean(cat(3,tinfo.delay),3);
                            end
                        end
                        
                    case {'central'}
                        C=cat(1,output(dd).condition(kk).data{:});
                        for ww=1:length(tinx{1,kk}) %for each window
                            windows(ww).mean=temp(ww).est;
                            windows(ww).median=tinfo(ww).median;
                            windows(ww).mode=tinfo(ww).mode;
                            windows(ww).std=tinfo(ww).std;
                            windows(ww).samples=tinfo(ww).samples;
                            windows(ww).missingdata=tinfo(ww).missingdata;
                        end
                        stemp=struct('label',tinx(2,kk),...
                            'mean', nanmean(C), ...
                            'median', nanmedian(C),...
                            'mode', mode(C),...
                            'std', nanstd(C),...
                            'samples', sum([tinfo.samples]), ...
                            'missingdata', sum(cat(1,tinfo.missingdata),1),...
                            'windows', windows);
                        
                    case {'auc'}
                        w=cat(1,tinfo.samples)-cat(1,tinfo.missingdata); %samples
                        C=cat(1,temp.est); %mean per window
                        for ww=1:length(tinx{1,kk}) %for each window
                            windows(ww).est=temp(ww).est;
                            windows(ww).samples=tinfo(ww).samples;
                            windows(ww).missingdata=tinfo(ww).missingdata;
                        end
                        stemp=struct('label',tinx(2,kk),...
                            'est', sum(C.*w,1)./sum(w,1), ...
                            'samples', sum([tinfo.samples]), ...
                            'missingdata', sum(cat(1,tinfo.missingdata),1),...
                            'windows', windows);
                        
                    case 'psd'
                        w=[];
                        w(1,1,:)=[tinfo.samples]; %samples
                        C=cat(3,temp.est); %psd per window
                        psd=sum(C.*w,3)./sum(w,3); %weighted mean
                        w=[tinfo.samples]'; %samples
                        C=cat(1,tinfo.totalpower);
                        tpower=sum(C.*w,1)./sum(w,1);
                        for ww=1:length(tinx{1,kk})
                            windows(ww).est=temp(ww).est;
                            windows(ww).samples=tinfo(ww).samples;
                            windows(ww).resampled=tinfo(ww).resampled;
                        end
                        stemp=struct('label',tinx(2,kk),...
                            'est', psd,...
                            'totalpower',tpower,...
                            'samples',sum([tinfo.samples]),...
                            'resampled',sum(cat(1,tinfo.resampled),1),...
                            'F',tinfo(1).F,...
                            'fullF',tinfo(1).fullF,...
                            'windows', windows);
                        
                    case 'dtw'
                        w=[];
                        w(1,1,:)=[tinfo.samples]; %samples
                        C=cat(3,temp.est);
                        for ww=1:length(tinx{1,kk})
                            windows(ww).est=temp(ww).est;
                            windows(ww).samples=tinfo(ww).samples;
                            windows(ww).missingdata=tinfo(ww).missingdata;
                        end
                        stemp=struct('label',tinx(2,kk),...
                            'est',sum(C.*w,3)./sum(w,3), ... %weighted mean
                            'samples', sum(w,3), ...
                            'missingdata', sum(cat(1,tinfo.missingdata),1), ...
                            'windows', windows);
                        
                    case 'scr'
                        temp=[temp.est];
                        w=cat(1,temp.CDAn); %weight by number of SCRs
                        C=cat(1,temp.CDAampu); %mean amplitude
                        CDAampu=sum(C.*w,1)./sum(w,1); CDAampu(isnan(CDAampu))=0;
                        C=cat(1,temp.CDAampstd); %std amplitude
                        CDAampstd=sum(C.*w,1)./sum(w,1); CDAampstd(isnan(CDAampstd))=0;
                        w=cat(1,temp.TTPn); %weight by number of TTPs
                        C=cat(1,temp.TTPampu); %mean amplitude
                        TTPampu=sum(C.*w,1)./sum(w,1); TTPampu(isnan(TTPampu))=0;
                        C=cat(1,temp.TTPampstd); %std amplitude
                        TTPampstd=sum(C.*w,1)./sum(w,1); TTPampstd(isnan(TTPampstd))=0;
                        w=cat(1,temp.samples)-cat(1,temp.missingdata);
                        C=cat(1,temp.CDAppm);
                        CDAppm=sum(C.*w,1)./sum(w,1);
                        C=cat(1,temp.TTPppm);
                        TTPppm=sum(C.*w,1)./sum(w,1);
                        if sum([tinfo.xwin])==0 %no windowing applied
                            CDAtseries=[];
                            TTPtseries=[];
                            timeline=[];
                        else %windowed
                            timeline=downsample(time,sr*tinfo(1).xstep);
                            CDAtseries=nan(length(timeline),subj);
                            TTPtseries=nan(length(timeline),subj);
                            for ww=1:length(tinx{1,kk}) %for each window
                                [~,t1]=min(abs(timeline-tinx{1,kk}{ww}(1))); %start idx
                                CDAtseries(t1:t1+length(temp(ww).CDAppmt(:,1))-1,:) = ...
                                    temp(ww).CDAppmt;
                                TTPtseries(t1:t1+length(temp(ww).TTPppmt(:,1))-1,:) = ...
                                    temp(ww).TTPppmt;
                            end
                        end
                        for ww=1:length(tinx{1,kk})
                            temp(ww).inx=tinx{1,kk}{ww};
                            temp(ww).t=tinx{1,kk}{ww};
                        end
                        stemp=struct('label',tinx(2,kk),...
                            'CDAn',sum(cat(1,temp.CDAn),1), ...
                            'CDAampu',CDAampu,...
                            'CDAampstd',CDAampstd,...
                            'CDAppm',CDAppm,...
                            'CDAtseries',CDAtseries,...
                            'TTPn',sum(cat(1,temp.TTPn),1), ...
                            'TTPampu',TTPampu,...
                            'TTPampstd',TTPampstd,...
                            'TTPppm',TTPppm,...
                            'TTPtseries',TTPtseries,...
                            'timeline',timeline,...
                            'samples', sum([temp.samples]), ...
                            'missingdata', sum(cat(1,temp.missingdata),1), ...
                            'windows', temp);
                end
                if kk==1
                    output(dd).analysis(aa).condition=stemp;
                else
                    output(dd).analysis(aa).condition= [ output(dd).analysis(aa).condition ...
                        stemp];
                end
%             catch
%                 zargs=cellfun(@num2str,args,'un',0);
%                 warning(['Major error in analysis [Arguments:' zargs{:} ]);
% 
%             end
        end
    end
end
fname=strsplit(pewtr3.file.sdir,{'\','/'});
if p==1, delete(f); end %close waitbar
pewtr4=struct('session',fname{end}, 'output', output, 'time', pewtr3.time, ...
    'subj', {pewtr3.subj}, 'file', pewtr3.file, 'groups', {pewtr3.bgfs}, ...
    'conditions', pewtr3.conditions, 'covar', {pewtr3.covar});
disp('Analysis complete')

