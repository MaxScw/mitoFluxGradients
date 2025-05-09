function DownsampledDecay_struct=DownsampleDecay(decay_struct,NumPhotonsOutput,Channel)
%DownsampleDecay takes an input decay struct and returns a new decay struct
%where a random subset of the photons have been chosen to create a new,
%downsampled decay.
%Inputs:    decay_struct: the source decay_struct file
%           NumPhotonsOutput: how many photons should be in the new decay
%           Channel: Which channel to use. 1-joint; 2-mito; 3-cyto
%Outputs:   DownsamplesDecay_struct: the downsampled decay_struct file
%           which only has NumPhotonsOutput photons in it, all the rest is the same
%Written by Will Conway, Needleman Lab, January 23, 2020


%Extract the decay you want
ExtractedDecay=decay_struct.decay(:,Channel);

%Generate a photon list
PhotonList=[];
for Index=1:length(ExtractedDecay)
    if ExtractedDecay(Index)==0
    else
        PhotonList=[PhotonList; Index*ones(ExtractedDecay(Index),1)];
    end
end

%Decide which photons to pull
Indexes=randperm(length(PhotonList));
DownsampledPhotonList=PhotonList(Indexes(1:NumPhotonsOutput));

%Rebin all of the photons
DownsampledDecay=hist(DownsampledPhotonList,1:1:256);

%Save the downsampled decay back into the decay_struct 
DownsampledDecay_struct=decay_struct;
DownsampledDecay_struct.decay(:,Channel)=DownsampledDecay;

end