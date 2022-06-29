function [output,varargout] = p_filtering(data,varargin)
%% Data filtering function for PEWTR
% output = pfiltering(data,'FilterType',...)
% data: data input is a matrix of data to be filtered, where each column
% represents a single continuous signal to be filtered.
% Replace 'FilterType' with filter of choice and necessary following
% arguments.
%
% pfiltering(data,'SG',Order,Window,Samplerate) 
% Savitzky-Golay filter. 
% Order: numeric value of polynomial order for filter.
% Window: numeric value of length of filter in seconds.
% Samplerate: numeric value for sample rate of signal in Hz.
%
% pfiltering(data,'Butter',Order,Frequency,Samplerate)
% Lowpass digital butterworth filter.
% Order: numeric value of order for filter.
% Frequency: numeric value of cutoff frequency.
% Samplerate: numeric value for sample rate of signal in Hz.

if numel(varargin)<1 %if no input arguments
    return
end

if ~isempty(find(strcmpi(varargin,'SG'), 1))
    %%S-G filter
    inx=find(strcmpi(varargin,'SG'), 1);
    win=varargin{inx+2}*varargin{inx+3};
    check=0;
    if rem(win,2)==0
        win=win+1;
        check=1;
    end
    
    if rem(win,1)~=0
        win=ceil(win);
        check=1;
    end
    
    if check==1
        errordlg(['Savitzky-Golay filtering requires odd number of samples '...
            'for window length. Window adjusted to ' num2str(win) ' samples ('...
            num2str(win/varargin{inx+3}) ' seconds)'],'Warning - Window length adjustment');
    end

    output = sgolayfilt(data,varargin{inx+1},win);
    varargout{1}=varargin{inx+1};
    varargout{2}=win/varargin{inx+3};
end

if ~isempty(find(strcmpi(varargin,'butter'), 1))
   %% Lowpass Butterworth filter
   inx=find(strcmpi(varargin,'butter'), 1);
   [b,a] = butter(varargin{inx+1},varargin{inx+2}/(varargin{inx+3}/2));
   output=filter(b,a,data);
   varargout{1}=varargin{inx+1};
   varargout{2}=varargin{inx+2};
    
end
