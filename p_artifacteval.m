function [varargout]=p_artifacteval(data,hz,varargin)

%% Evaluation metrics for EDA removal settings

% jump(jumpwin,jumpchange,slopeSD); short(win); entropy(win)

[row,col]=size(data);
varargout={};

%% Counts of Jump Artifacts (C1)
%Absolute slope between an EDA datapoint and the following four datapoints
%is computed then averaged. Threshold should be set as the minimum between
%the value of 2SD of the EDA signal and the equivalent to the change of
%10uS per second

if sum(strcmpi(varargin,'jump'))~=0
    inx=find(strcmpi(varargin,'jump'));
    
    ds=hz*(varargin{inx+1}/1000);
    if ceil(ds)~=ds
        errordlg(['Selected window length for slope calculation is too large for' ...
            ' sampling rate. Window length adjusted to ' num2str(ceil(ds)/hz) ...
            'ms (' num2str(ceil(ds)) ' data points)'],'Artifact Eval Warning');
        ds=ceil(ds);
        varargin{inx+1}=(ds/hz)*1000;
    end
    jcount=zeros(size(data));
    jcount(end-ds+1:end,:)=[];
    for ii=1:col %for each subject
        for jj=1:numel(jcount(:,ii)) %for each datapoint
            jcount(jj,ii)=mean(abs(data(jj,ii)-data(jj+1:jj+ds,ii)),'omitnan');
        end
    end
    
    varargout{end+1}=sum(jcount>((varargin{inx+1}/1000)*varargin{inx+2}));
    
    % jmean=mean(jcount,'omitnan');
    M1=NaN(1,col);
    M2=NaN(1,col);
    jsd=std(jcount,'omitnan');
    for ii=1:col
        jinx=jcount(:,ii)>=(jsd(ii)*varargin{inx+3});
        M1(ii)=mean(jcount(jinx,ii));
        M2(ii)=sum(jcount(jinx,ii));
    end
    
    varargout{end+1}=M1;
    varargout{end+1}=M2;
end

%% Counts of Short Epochs
% count the number of remaining EDA segments shorter than 10s

if sum(strcmpi(varargin,'short'))~=0
    inx=find(strcmpi(varargin,'short'));
    
    C2=zeros(1,col);
    M2=zeros(1,col);
    y=diff(~isnan(vertcat(NaN(1,col),data,NaN(1,col)))); %pad, ~isnan, get diff
    for ii=1:col %for each subject
        x=find(y(:,ii)==-1)-find(y(:,ii)==1); %get length of each continuous epoch
        C2(ii)=sum(x<(varargin{inx+1}*hz)); %Count epochs less than specified length
        M2(ii)=sum(x(x<(varargin{inx+1}*hz)))/numel(~isnan(data(:,ii))); %Percentage of total non-zero data
    end
    varargout{end+1}=C2;
    varargout{end+1}=M2;
end

%% Decreased Entropy
% In each 10 second window, Shannon entropy based on probability
% distribution using the histogram technique (average entropy per
% participant)

if sum(strcmpi(varargin,'entropy'))~=0
    inx=find(strcmpi(varargin,'entropy'));

y=normalize(data);
M2=NaN(row,col);
for ii=1:col
    for jj=1:varargin{inx+1}*hz:(row-rem(row,varargin{inx+1}*hz))
        M2(jj,ii)=wentropy(y(jj:jj+varargin{inx+1}*hz,ii),'shannon');
    end
end
varargout{end+1}=mean(M2,'omitnan');
end
 
