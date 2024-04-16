function [scr,info] = p_scr(SCR,win,varargin)
% Extracts SCR measures from pewtr3.SCR structure
%
% [scr,info]=p_scr(SCR,win,meanwin,meanstep)
%
% INPUTS:
% SCR - pewtr3.SCR structure containing CDA and TTP fields
%
% WIN - start and end times (in minutes) for time window in sochro timeline 
% to extract SCR measures e.g. [0 10] for 0 to 10 minutes
%
% OPTIONAL INPUTS:
%
% MEANWIN - calculate SCR frequency using a moving window specified (in 
% seconds) by MEANWIN. If unspecified, mean SCR rate per minute is used.
%
% MEANSTEP - step size (in seconds) used for calculating SCR rate per
% window

if isempty(varargin)
    xwin=0;
    xstep=0;
else
    xwin=varargin{1};
    xstep=varargin{2};
end
info.win=win; info.xwin=xwin; info.xstep=xstep;
xwin=xwin/60; xstep=xstep/60; %convert window and step to minutes
subj=length(SCR); %replaced length([SCR.CDA])

scr=struct('CDAn',zeros(1,subj),'CDAampu',zeros(1,subj),'CDAampstd',...
    zeros(1,subj),'CDAppm',zeros(1,subj),'CDAtime',...
    win(1):xstep:win(2),'CDAppmt',zeros(length(win(1):xstep:win(2)),subj),...
    'TTPn',zeros(1,subj),'TTPampu',zeros(1,subj),'TTPampstd',zeros(1,subj),...
    'TTPppm',zeros(1,subj),'TTPtime',win(1):xstep:win(2),...
    'TTPppmt',zeros(length(win(1):xstep:win(2)),subj));

for ii=1:subj
    %% CDA
    if isempty(SCR(ii).CDA)
        scr.CDAn(ii)=NaN;
        scr.CDAampu(ii)=NaN;
        scr.CDAampstd(ii)=NaN;
        scr.CDAppm(ii)=NaN;
        scr.CDAppmt(ii)=NaN;
    else
        inx=find(SCR(ii).CDA.onset>=win(1) & SCR(ii).CDA.onset<=win(2));
        if ~isempty(inx)
            scr.CDAn(ii)=length(inx);
            scr.CDAampu(ii)=mean(SCR(ii).CDA.amp(inx));
            scr.CDAampstd(ii)=std(SCR(ii).CDA.amp(inx));
            scr.CDAppm(ii)=length(inx)/(win(2)-win(1));
            if xwin~=0
                for tt=1:length(scr.CDAtime)
                    tinx=find(SCR(ii).CDA.onset>=scr.CDAtime(tt) ...
                        & SCR(ii).CDA.onset<=(scr.CDAtime(tt)+xwin));
                    if isempty(tinx)
                        scr.CDAppmt(tt,ii)=0;
                    else
                        scr.CDAppmt(tt,ii)=length(tinx)/xwin;
                    end
                end
            end
        end
    end
    
    %% TTP
    if isempty(SCR(ii).TTP)
        scr.TTPn(ii)=NaN;
        scr.TTPampu(ii)=NaN;
        scr.TTPampstd(ii)=NaN;
        scr.TTPppm(ii)=NaN;
        scr.TTPppmt(ii)=NaN;
    else
        inx=find(SCR(ii).TTP.onset>=win(1) & SCR(ii).TTP.onset<=win(2));
        if ~isempty(inx)
            scr.TTPn(ii)=length(inx);
            scr.TTPampu(ii)=mean(SCR(ii).TTP.amp(inx));
            scr.TTPampstd(ii)=std(SCR(ii).TTP.amp(inx));
            scr.TTPppm(ii)=length(inx)/(win(2)-win(1));
            if xwin~=0
                for tt=1:length(scr.TTPtime)
                    tinx=find(SCR(ii).TTP.onset>=scr.TTPtime(tt) ...
                        & SCR(ii).TTP.onset<=(scr.TTPtime(tt)+xwin));
                    if isempty(tinx)
                        scr.TTPppmt(tt,ii)=0;
                    else
                        scr.TTPppmt(tt,ii)=length(tinx)/xwin;
                    end
                end
            end
        end
    end
end