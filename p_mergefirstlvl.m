function [pewtr] = p_mergefirstlvl(pewtr4)
% Merges analysis substructures according to session label
% Input should be a concatonated pewtr4 structure

session=unique({pewtr4.session}); %get unique session labels
for ss=1:length(session)
    inx=find(strcmpi(session(ss),{pewtr4.session})); %indices of session ss
    pewtr(ss)=pewtr4(inx(1)); %copy info to first index
    pewtr(ss).output=[];
    output=[pewtr4(inx).output]; %concatonate analysis outputs
    datatype=unique({output.datatype}); %indices of unique datatypes
    for dd=1:length(datatype) %for each datatype
        pewtr(ss).output(dd).datatype=datatype{dd};
        inx=find(strcmpi(datatype(dd),{output.datatype}));
        pewtr(ss).output(dd).analysis=[output(inx).analysis];
        pewtr(ss).output(dd).condition=[output(inx).condition];
        %remove empty analysis slots
        dinx=[];
        for aa=1:length(pewtr(ss).output(dd).analysis)
            if isempty(pewtr(ss).output(dd).analysis(aa).label)
                dinx=[dinx aa];
            end
        end
        pewtr(ss).output(dd).analysis(dinx)=[];
        %remove duplicate conditions
        [~,inx]=unique({pewtr(ss).output(dd).condition.label});
        pewtr(ss).output(dd).condition=pewtr(ss).output(dd).condition(inx);    
    end
end

