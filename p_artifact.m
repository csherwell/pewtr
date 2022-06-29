function [output,varargout] = p_artifact(data,hz,varargin)

% output = p_artifact(_____,'temp',tempdata,temphz,[min max])
% output = p_artifact(_____,'range',[min max])
% output = p_artifact(_____,'jump',jumpwin,jumpchange)
% output = p_artifact(_____,'trans',transwin)
% output = p_artifact(_____,'interloss',min)


% p_artifact(data,4,'temp',dispdat.TEMP,4,[30 40],'range',[0.05 60],'jump',250,10,'trans',1,'interloss',5)
%% Artifact detection and evaluation metrics for EDA

output=zeros([size(data),5]); %output is logical array of identified artifacts
[~,col]=size(data); %number of participants in data

%% 1. Temperature in range
if sum(strcmpi('temp',varargin))~=0
   inx=find(strcmpi('temp',varargin));
   if hz~=varargin{inx+2}
       errordlg(['Warning - PEWTR cannot run temperature range exclusion '...
           'where EDA and temperature signals are recorded at different '...
           'sampling frequencies. Artifact detection aborted'],'Error');
       return
   end

   d1=varargin{inx+1}<varargin{inx+3}(1); %min
   d2=varargin{inx+1}>varargin{inx+3}(2); %max
   d3=isnan(varargin{inx+1}); %nan
   output(:,:,1)=logical(d1+d2+d3);
end

%% 2. EDA range
if sum(strcmpi('range',varargin))~=0
   inx=find(strcmpi('range',varargin));
   d1=data<varargin{inx+1}(1); %min
   d2=data>varargin{inx+1}(2); %max
   output(:,:,2)=logical(d1+d2);
end

%% 3. Jump artifacts
if sum(strcmpi('jump',varargin))~=0
    inx=find(strcmpi('jump',varargin));
    ds=hz*(varargin{inx+1}/1000); %window size in data samples calculation
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
    jcount=vertcat(jcount,zeros(ds,col));
    
    d1=jcount>(varargin{inx+1}/1000)*varargin{inx+2};
    d1=vertcat(d1,NaN(ds-1,col));
    output(:,:,3)=logical(d1);
end

%% 4. Transitional Period
if sum(strcmpi('trans',varargin))~=0
   inx=find(strcmpi('trans',varargin));
   ds=hz*varargin{inx+1}; %get transitional period window in samples
   P=sum(output,3)==0; %Combine all rejected data, 1=present
   P(isnan(data))=0; %Remove missing data
   
   d1=ds>=movsum(P,[ds 0]);
   d2=ds>=movsum(P,[0 ds]);
   out=logical(d1+d2);
   out(~P)=0;
   output(:,:,4)=out;
end

%% 5. Inter-loss Period

if sum(strcmpi('interloss',varargin))~=0
   inx=find(strcmpi('interloss',varargin));
   ds=hz*varargin{inx+1};
   P=sum(output,3)==0; %Combine all rejected data, 1=present
   P(isnan(data))=0; %Remove missing data

   y=diff(vertcat( zeros(1,col), P, zeros(1,col))); %pad, ~isnan, diff
   for ii=1:col
       x1=find(y(:,ii)==1); %start
       x2=find(y(:,ii)==-1); %end
       x3=x2-x1; %epoch length
       z=find(x3<ds); %indices of subthreshold epochs
       if ~isempty(z)
          for jj=1:numel(z) %for each epoch to reject
             output(x1(z(jj)):x2(z(jj))-1,ii,5)=1; 
          end
       end
   end
end


%% Evaluation of missing data

M=sum(isnan(data)); %get sum of missing data
D=sum(~isnan(data)); %get sum of existing data
A=sum((sum(output,3)~=0)-(isnan(data))); %sum of artifact (minus already missing data)
varargout={D,M,A};

