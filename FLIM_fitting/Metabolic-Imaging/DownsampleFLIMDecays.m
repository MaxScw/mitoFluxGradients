function DownsampleFLIMDecays(acqpath)
%DownsampleFLIMDecays loads all of the decay outputs from
%FLIM_decay_from_probmap and downsamples the decays so that the ICM and the
%Troph have equal numbers of photons.
%Inputs:    -acqpath: The path where all of the files are saved
%Outputs:   -Overwrites all of the decay files with the new downsampled
%           files
%           -Saves all of the old decays in a folder called OriginalDecays
%           if there are no previous original decays
%Written by Will Conway, Needleman Lab January 24, 2020

cd(acqpath)
if exist('OriginalDecays', 'dir')
    cd OriginalDecays
    PreviousSavedDecays=1;
else
    PreviousSavedDecays=0;
end
Listing=dir(pwd);

%Load all of the decays in

LoadIndex=1;
for ListingIndex=3:length(Listing)
    Filename=Listing(ListingIndex).name;
    if Filename(1:6)=='decays'
        LoadedDecayNames{LoadIndex}=Filename;
        LoadedDecay(LoadIndex)=load(fullfile(acqpath,Filename));
        LoadIndex=LoadIndex+1;
    else
    end
    
end



%Save a backup
if PreviousSavedDecays==0
    mkdir('OriginalDecays')
    cd OriginalDecays
    
    for FileIndex=1:length(LoadedDecay)
        SaveName=LoadedDecayNames{FileIndex};
        SaveFile=LoadedDecay(FileIndex);
        decay_struct=SaveFile.decay_struct;
        save(SaveName,'decay_struct')
    end
   
else
end

cd ..

%Downsample all of the decays
Channel=1;
for PositionIndex=0:round(length(LoadedDecay))-1
    
    %Pull the decays for the position
    PositionDecays=LoadedDecay(~cellfun('isempty',regexp(LoadedDecayNames,strcat('Pos',num2str(PositionIndex)))));
    PositionDecayNames={LoadedDecayNames{~cellfun('isempty',regexp(LoadedDecayNames,strcat('Pos',num2str(PositionIndex))))}};
    
    %Pull the ICM and TrophDecays
    TrophDecay=PositionDecays(~cellfun('isempty',regexp(PositionDecayNames,'Troph')));
    if ~isempty(TrophDecay)
            TrophDecay=TrophDecay.decay_struct;
    else
    end
    
    ICMDecay=PositionDecays(~cellfun('isempty',regexp(PositionDecayNames,'ICM')));
    if ~isempty(ICMDecay)
    ICMDecay=ICMDecay.decay_struct;
    else
    end
    
     if ~isempty(ICMDecay)
          if ~isempty(TrophDecay)
            for DecayIndex=1:length(TrophDecay)

                if isempty(TrophDecay{DecayIndex})
                    
                else
                     if isempty(ICMDecay{DecayIndex})
                     else
                    %Compute how many photons to downsample
                    NumberPhotonsTroph=sum(TrophDecay{DecayIndex}.decay(:,Channel));
                    NumberPhotonsICM=sum(ICMDecay{DecayIndex}.decay(:,Channel));
                    NumberPhotonsToSample=min(NumberPhotonsTroph,NumberPhotonsICM);

                    %Downsample the decays
                    TrophDecay{DecayIndex}=DownsampleDecay(TrophDecay{DecayIndex},NumberPhotonsToSample,Channel);
                    ICMDecay{DecayIndex}=DownsampleDecay(ICMDecay{DecayIndex},NumberPhotonsToSample,Channel);
                     end
                end

            end
          else
          end
     else
     end
    
    %Put the new ICM and Troph Decays back in the original positions
    TempDecay.decay_struct=TrophDecay;
    PositionDecays(~cellfun('isempty',regexp(PositionDecayNames,'Troph')))=TempDecay;
    TempDecay2.decay_struct=ICMDecay;
    PositionDecays(~cellfun('isempty',regexp(PositionDecayNames,'ICM')))=TempDecay2;
    
    LoadedDecay(~cellfun('isempty',regexp(LoadedDecayNames,strcat('Pos',num2str(PositionIndex)))))=PositionDecays;
    
end

%Save a downsampled decays
clear decay_struct
for FileIndex=1:length(LoadedDecay)
   
    SaveName=LoadedDecayNames{FileIndex};
    SaveFile=LoadedDecay(FileIndex);
    decay_struct=SaveFile.decay_struct;
    save(SaveName,'decay_struct')
    
end




 end