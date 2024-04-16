function [pewtr2] = p_importempatica(sdir)
%% PEWTR - import empatica data
% INPUT: sdir = directory containing directories of participant data
% OUTPUT: PEWTR formatted data structure

subdir=dir(sdir); %get subdirectories in main directory
subdir={subdir([subdir.isdir]==1).name}; %get names (only of directories)
d=strncmp('.',subdir,1); %get indices of non-directories
subdir(d)=[]; %remove them

%% Collect individual subject data

fn=zeros(1,numel(subdir)); %index for subjects to remove
urdata=struct([]);
for ii = 1:numel(subdir) %for each directory    
    urdata(ii).urEDA=[]; urdata(ii).urACC=[];
    urdata(ii).urBVP=[]; urdata(ii).urIBI=[];
    urdata(ii).urHR=[]; urdata(ii).urTEMP=[];
    
    fl=dir([sdir filesep subdir{ii}]); %get file list of everything in directory
    disp(['Loading subject folder: ' subdir{ii}])
    
    for kk=1:numel(fl) %for every file in the file list of subdirectory ii
        % LOAD EDA
        if strcmp('EDA.csv',fl(kk).name)%if file kk is an EDA.csv file
            urdata(ii).urEDA=importdata([sdir filesep subdir{ii} filesep 'EDA.csv']); %import raw data
            if length(urdata(ii).urEDA)<3 %if there is only header info
                urdata(ii).urEDA=[];
            end
            
        % LOAD ACC
        elseif strcmp('ACC.csv',fl(kk).name)%if file kk is an ACC.csv
            urdata(ii).urACC=importdata([sdir filesep subdir{ii} filesep 'ACC.csv']);
            
        % LOAD BVP
        elseif strcmp('BVP.csv',fl(kk).name)%if file kk is a BVP.csv
            urdata(ii).urBVP=importdata([sdir filesep subdir{ii} filesep 'BVP.csv']);
            if length(urdata(ii).urBVP)<3 %if there is only header info
                urdata(ii).urBVP=[];
            end
            
        % LOAD IBI
        elseif strcmp('IBI.csv',fl(kk).name)%if file kk is a IBI.csv
            temp=importdata([sdir filesep subdir{ii} filesep 'IBI.csv']);
            if ~isempty(temp)
                try
                    urdata(ii).urIBI(1,1)=str2double(temp.textdata{1,1});
                    urdata(ii).urIBI(1,2)=0;
                    urdata(ii).urIBI(2:length(temp.data(:,1))+1,1)=temp.data(:,1)+urdata(ii).urIBI(1,1);
                    urdata(ii).urIBI(2:length(temp.data(:,1))+1,2)=temp.data(:,2);
                catch
                end
            end
            
            % LOAD HR
        elseif strcmp('HR.csv',fl(kk).name)%if file kk is a HR.csv
            urdata(ii).urHR=importdata([sdir filesep subdir{ii} filesep 'HR.csv']);
            if length(urdata(ii).urHR)<3 %if there is only header info
                urdata(ii).urHR=[];
            end
            
            % LOAD TEMP
        elseif strcmp('TEMP.csv',fl(kk).name) %if file kk is a TEMP.csv
            urdata(ii).urTEMP=importdata([sdir filesep subdir{ii} filesep 'TEMP.csv']);
            if length(urdata(ii).urTEMP)<3 %if there is only header info
                urdata(ii).urTEMP=[];
            end
        end
        
    end %end of file loop kk

    if isempty(urdata(ii).urEDA) && isempty(urdata(ii).urACC) && ...
            isempty(urdata(ii).urBVP) && isempty(urdata(ii).urIBI) && ...
            isempty(urdata(ii).urHR) && isempty(urdata(ii).urTEMP) 
        disp(['NO DATA FOUND IN DIR ' subdir{ii} ...
            ' - removing from participant list']);
        fn(ii)=1;
    end
    
end

urdata(logical(fn)) = []; %Remove missing datasets
subj=subdir(~(logical(fn)));

disp('Formatting data');

%% Get global timing parameters
globalon=nan(length(urdata),6); 
globaloff=nan(length(urdata),6);
IBIc=0;
for ii=1:length(urdata)
    %EDA
    if ~isempty(urdata(ii).urEDA)
        urdata(ii).EDAhz=urdata(ii).urEDA(2,1); %log sampling rate
        urdata(ii).EDAon=urdata(ii).urEDA(1,1); %log EDA start time
        urdata(ii).EDAoff= urdata(ii).EDAon + ...
            (length(urdata(ii).urEDA(3:end,1))/urdata(ii).EDAhz); %log EDA end time as start + (length / samplerate)
        globalon(ii,1)=urdata(ii).EDAon;
        globaloff(ii,1)=urdata(ii).EDAoff;
    else
        urdata(ii).EDAhz=[];
    end
    
    %BVP
    if ~isempty(urdata(ii).urBVP)
        urdata(ii).BVPhz=urdata(ii).urBVP(2,1); %log sampling rate
        urdata(ii).BVPon=urdata(ii).urBVP(1,1); %log BVP start time
        urdata(ii).BVPoff= urdata(ii).BVPon + ...
            (length(urdata(ii).urBVP(3:end,1))/urdata(ii).BVPhz); %log BVP end time as start + (length / samplerate)
        globalon(ii,2)=urdata(ii).BVPon;
        globaloff(ii,2)=urdata(ii).BVPoff;
    else
        urdata(ii).BVPhz=[];
    end
    
    %ACC
    if ~isempty(urdata(ii).urACC)
        urdata(ii).ACChz=urdata(ii).urACC(2,1); %log sampling rate
        urdata(ii).ACCon=urdata(ii).urACC(1,1); %log ACC start time
        urdata(ii).ACCoff= urdata(ii).ACCon + ...
            (length(urdata(ii).urACC(3:end,1))/urdata(ii).ACChz); %log ACC end time as start + (length / samplerate)
        globalon(ii,3)=urdata(ii).ACCon;
        globaloff(ii,3)=urdata(ii).ACCoff;
    else
        urdata(ii).ACChz=[];
    end
    
    %IBI
    if ~isempty(urdata(ii).urIBI)
        urdata(ii).IBIon=urdata(ii).urIBI(1,1); %log IBI start time
        urdata(ii).IBIoff=urdata(ii).urIBI(end,1); %log IBI end time
        globalon(ii,4)=urdata(ii).IBIon;
        globaloff(ii,4)=urdata(ii).IBIoff;
        IBIc=1;
    end
    
    %HR
    if ~isempty(urdata(ii).urHR)
        urdata(ii).HRhz=urdata(ii).urHR(2,1); %log sampling rate
        urdata(ii).HRon=urdata(ii).urHR(1,1); %log BVP start time
        urdata(ii).HRoff= urdata(ii).HRon + ...
            (length(urdata(ii).urHR(3:end,1))/urdata(ii).HRhz); %log BVP end time as start + (length / samplerate)
        globalon(ii,5)=urdata(ii).HRon;
        globaloff(ii,5)=urdata(ii).HRoff;
    else
        urdata(ii).HRhz=[];
    end
    
    %TEMP
    if ~isempty(urdata(ii).urTEMP)
        urdata(ii).TEMPhz=urdata(ii).urTEMP(2,1); %log sampling rate
        urdata(ii).TEMPon=urdata(ii).urTEMP(1,1); %log TEMP start time
        urdata(ii).TEMPoff= urdata(ii).TEMPon + ...
            (length(urdata(ii).urTEMP(3:end,1))/urdata(ii).TEMPhz); %log TEMP end time as start + (length / samplerate)
        globalon(ii,6)=urdata(ii).TEMPon;
        globaloff(ii,6)=urdata(ii).TEMPoff;
    else
        urdata(ii).TEMPhz=[];
    end
end

%Find global times
globalon= min(min(globalon));
globaloff= max(max(globaloff));

%% EDA processing
if ~isempty([urdata.EDAhz])
    %Samplerate consistency check
    if std([urdata.EDAhz])~=0 %if there is variance in samplerates
        ButtonName = questdlg(['Inconsistent sample rates between EDA files. ' ...
            'Do you wish to resample at the lowest sampling rate of ',...
            num2str(min([urdata.EDAhz])) 'Hz?','Warning','Yes','No','Yes']);
        switch ButtonName
            case 'Yes'
                sr=min([urdata.EDAhz]);
                for ii=1:length(urdata)
                    if urdata(ii).EDAhz>sr
                        urdata(ii).urEDA=resample(urdata(ii).urEDA(3:end),sr,urdata(ii).EDAhz);
                    end
                end
            case 'No'
                file.Hz.EDAhz=-1;
                errordlg(['Inconsistent sample rates between EDA files. ' ...
                    'Errors in computations and data display may occur' ]);
        end
    else
        file.Hz.EDAhz=min([urdata.EDAhz]); %copy samplerate to global parameters
    end
    
    %Temporally align subjects biometrics and pad with NaN for equal vector length
    timeline= globalon:1/max([urdata.EDAhz]):globaloff; %calculate timeline from minimum and maximum stamps
    EDA=NaN(length(timeline),length(urdata)); %create empty sample*subject matrix
    for ii=1:length(urdata)
        if ~isempty(urdata(ii).urEDA)
            [~,inx]=min(abs(timeline-urdata(ii).EDAon)); %find closest matching starting point in timeline for subject ii
            EDA(inx:inx-1+length(urdata(ii).urEDA(3:end,1)),ii)= ... %assign data according to starting point in timeline
                urdata(ii).urEDA(3:end,1);
        end
    end
    time.EDA=timeline;
    time.EDAmins=((0:1/max([urdata.EDAhz]):timeline(end)-timeline(1)))./60;
end

%% BVP processing
if ~isempty([urdata.BVPhz])
    %Samplerate consistency check
    if std([urdata.BVPhz])~=0 %if there is variance in samplerates
        ButtonName = questdlg(['Inconsistent sample rates between BVP files. ' ...
            'Do you wish to resample at the lowest sampling rate of ',...
            num2str(min([urdata.BVPhz])) 'Hz?','Warning','Yes','No','Yes']);
        switch ButtonName
            case 'Yes'
                sr=min([urdata.BVPhz]);
                for ii=1:length(urdata)
                    if urdata(ii).BVPhz>sr
                        urdata(ii).urBVP=resample(urdata(ii).urBVP(3:end),sr,urdata(ii).BVPhz);
                    end
                end
            case 'No'
                file.Hz.BVPhz=-1;
                errordlg(['Inconsistent sample rates between BVP files. ' ...
                    'Errors in computations and data display may occur' ]);
        end
    else
        file.Hz.BVPhz=min([urdata.BVPhz]); %copy samplerate to global parameters
    end
    
    %Temporally align subjects biometrics and pad with NaN for equal vector length
    timeline= globalon:1/max([urdata.BVPhz]):globaloff; %calculate timeline from minimum and maximum stamps
    BVP=NaN(length(timeline),length(urdata)); %create empty sample*subject matrix
    for ii=1:length(urdata)
        if ~isempty(urdata(ii).urBVP)
        [~,inx]=min(abs(timeline-urdata(ii).BVPon)); %find closest matching starting point in timeline for subject ii
        BVP(inx:inx-1+length(urdata(ii).urBVP(3:end,1)),ii)= ... %assign data according to starting point in timeline
            urdata(ii).urBVP(3:end,1);
        end
    end
    time.BVP=timeline;
    time.BVPmins=((0:1/max([urdata.BVPhz]):timeline(end)-timeline(1)))./60;
end

%% IBI processing
if IBIc==1
    %Temporally align subjects biometrics and pad with NaN for equal vector length
    if isempty(max([urdata.BVPhz]))
        file.Hz.IBIhz=64;
        warndlg('No blood pulse volume data found: default inter-beat interval sampling rate is set to 64Hz');
    else
        file.Hz.IBIhz=max([urdata.BVPhz]);
    end
    
    timeline= globalon:1/file.Hz.IBIhz:globaloff; %calculate timeline from minimum and maximum stamps
    IBI=NaN(length(timeline),length(urdata)); %create empty sample*subject matrix
    for ii=1:length(urdata)
        if ~isempty(urdata(ii).urIBI)
            for tt=1:length(urdata(ii).urIBI(2:end,1))
                [~,inx]=min(abs(timeline-urdata(ii).urIBI(tt+1,1))); %find closest matching timepoint in timeline for subject ii
                IBI(inx,ii)=urdata(ii).urIBI(tt+1,2);
            end
        end
    end
    time.IBI=timeline;
    time.IBImins=((0:1/file.Hz.IBIhz:timeline(end)-timeline(1)))./60;
end

%% ACC processing
if ~isempty([urdata.ACChz])
    %Samplerate consistency check
    if std([urdata.ACChz])~=0 %if there is variance in samplerates
        ButtonName = questdlg(['Inconsistent sample rates between ACC files. ' ...
            'Do you wish to resample at the lowest sampling rate of ',...
            num2str(min([urdata.ACChz])) 'Hz?','Warning','Yes','No','Yes']);
        switch ButtonName
            case 'Yes'
                sr=min([urdata.ACChz]);
                for ii=1:length(urdata)
                    if urdata(ii).ACChz>sr
                        urdata(ii).urACC=resample(urdata(ii).urACC(3:end,:),sr,urdata(ii).ACChz);
                    end
                end
            case 'No'
                file.Hz.ACChz=-1;
                errordlg(['Inconsistent sample rates between ACC files. ' ...
                    'Errors in computations and data display may occur' ]);
        end
    else
        file.Hz.ACChz=min([urdata.ACChz]); %copy samplerate to global parameters
    end
    
    %Temporally align subjects biometrics and pad with NaN for equal vector length
    timeline= globalon:1/max([urdata.ACChz]):globaloff; %calculate timeline from minimum and maximum stamps
    ACC=NaN(length(timeline),length(urdata)); %create empty sample*subject matrix
    for ii=1:length(urdata)
        if ~isempty(urdata(ii).urACC)
            
        %Convert to g
        g=(urdata(ii).urACC(3:end,:)/64);
        % ENMO calculation
        enmo=(sqrt(sum(g.*g,2)))-1;
        enmo(enmo<0)=0;
        
            
            
        [~,inx]=min(abs(timeline-urdata(ii).ACCon)); %find closest matching starting point in timeline for subject ii
%         ACC(inx:inx-1+length(urdata(ii).urACC(3:end,1)),ii)= ... %assign data according to starting point in timeline
%             sum(urdata(ii).urACC(3:end,:).^2,2);
         ACC(inx:inx-1+length(urdata(ii).urACC(3:end,1)),ii)=enmo; %assign data according to starting point in timeline

        end
    end
    time.ACC=timeline;
    time.ACCmins=((0:1/max([urdata.ACChz]):timeline(end)-timeline(1)))./60;
end

%% HR processing
if ~isempty([urdata.HRhz])
    %Samplerate consistency check
    if std([urdata.HRhz])~=0 %if there is variance in samplerates
        ButtonName = questdlg(['Inconsistent sample rates between HR files. ' ...
            'Do you wish to resample at the lowest sampling rate of ',...
            num2str(min([urdata.HRhz])) 'Hz?','Warning','Yes','No','Yes']);
        switch ButtonName
            case 'Yes'
                sr=min([urdata.HRhz]);
                for ii=1:length(urdata)
                    if urdata(ii).HRhz>sr
                        urdata(ii).urHR=resample(urdata(ii).urHR(3:end),sr,urdata(ii).HRhz);
                    end
                end
            case 'No'
                file.Hz.HRhz=-1;
                errordlg(['Inconsistent sample rates between HR files. ' ...
                    'Errors in computations and data display may occur' ]);
        end
    else
        file.Hz.HRhz=min([urdata.HRhz]); %copy samplerate to global parameters
    end
    
    %Temporally align subjects biometrics and pad with NaN for equal vector length
    timeline= globalon:1/max([urdata.HRhz]):globaloff; %calculate timeline from minimum and maximum stamps
    HR=NaN(length(timeline),length(urdata)); %create empty sample*subject matrix
    for ii=1:length(urdata)
        if ~isempty(urdata(ii).urHR)
        [~,inx]=min(abs(timeline-urdata(ii).HRon)); %find closest matching starting point in timeline for subject ii
        HR(inx:inx-1+length(urdata(ii).urHR(3:end,1)),ii)= ... %assign data according to starting point in timeline
            urdata(ii).urHR(3:end,1);
        end
    end
    time.HR=timeline;
    time.HRmins=((0:1/max([urdata.HRhz]):timeline(end)-timeline(1)))./60;

end

%% TEMP processing
if ~isempty([urdata.TEMPhz])
    %Samplerate consistency check
    if std([urdata.TEMPhz])~=0 %if there is variance in samplerates
        ButtonName = questdlg(['Inconsistent sample rates between TEMP files. ' ...
            'Do you wish to resample at the lowest sampling rate of ',...
            num2str(min([urdata.HRhz])) 'Hz?','Warning','Yes','No','Yes']);
        switch ButtonName
            case 'Yes'
                sr=min([urdata.TEMPhz]);
                for ii=1:length(urdata)
                    if urdata(ii).TEMPhz>sr
                        urdata(ii).urTEMP=resample(urdata(ii).urTEMP(3:end),sr,urdata(ii).TEMPhz);
                    end
                end
            case 'No'
                file.Hz.TEMPhz=-1;
                errordlg(['Inconsistent sample rates between TEMP files. ' ...
                    'Errors in computations and data display may occur' ]);
        end
    else
        file.Hz.TEMPhz=min([urdata.TEMPhz]); %copy samplerate to global parameters
    end
    
    %Temporally align subjects biometrics and pad with NaN for equal vector length
    timeline= globalon:1/max([urdata.TEMPhz]):globaloff; %calculate timeline from minimum and maximum stamps
    TEMP=NaN(length(timeline),length(urdata)); %create empty sample*subject matrix
    for ii=1:length(urdata)
        if ~isempty(urdata(ii).urTEMP)
        [~,inx]=min(abs(timeline-urdata(ii).TEMPon)); %find closest matching starting point in timeline for subject ii
        TEMP(inx:inx-1+length(urdata(ii).urTEMP(3:end,1)),ii)= ... %assign data according to starting point in timeline
            urdata(ii).urTEMP(3:end,1);
        end
    end
    time.TEMP=timeline;
    time.TEMPmins=((0:1/max([urdata.TEMPhz]):timeline(end)-timeline(1)))./60;
end

%% Introduce important file parameters
file.sdir=sdir;
file.unittime='UTC';
file.delay=0;

%fields
pewtr2.subj=subj;
pewtr2.file=file;
pewtr2.time=time;
pewtr2.conditions=cell(1);
pewtr2.covar=cell(1);
pewtr2.bgfs=cell(1);

%% Copy data if available
if exist('EDA','var'), pewtr2.data.EDA=EDA; end
if exist('BVP','var'), pewtr2.data.BVP=BVP; end
if exist('IBI','var'), pewtr2.data.IBI=IBI; end
if exist('ACC','var'), pewtr2.data.ACC=ACC; end
if exist('HR','var'), pewtr2.data.HR=HR; end
if exist('TEMP','var'), pewtr2.data.TEMP=TEMP; end

pewtr2.process=cell(numel(fieldnames(pewtr2.data)),10);

disp('Data import complete');

