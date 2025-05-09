
% 2014-08-28
% Changed to using 750 and 845nm, so change time is 16 instead of 25s. Use
% 22s wait time for laser. 
% Also using uMan DAC card to trigger and a single wait channel to absorb
% all the wait time.
%% Multi pos, z scan
z = 3;
clc
P = 8;
E = 20000; % This means acquistion time PLUS save time. Allow >5s save
if P>1 % include position change time if there are more than 1 positions
    lag = z*2*80 + 4*200 + 840 + 4000
    % lag for one position, dual scan.
    % 80ms per z change, 200ms to change filter wheel, 840ms to change
    % position
    % about 3 seconds to change the HWP
else
    lag = z*2*80 + 4*200 + 4000
end
% lagperwaitfr = lag/(z*2); % obsolete now with DAC
% waitperfr = 22000/z-lagperwaitfr*1.2; % obsolete now with DAC

lasertime = 23000+E*z;
uMantime = lasertime*2*P;
waitchannel = 23000-lag; % all lag time absorbed into one wait channel
[{'lag','lasertime','uMantime','waitchannel'};...
    {lag,lasertime,uMantime,waitchannel}]

%% dual timelapse, 1 pos
z = 1;
P = 2;
E = 5000;
if P>1 % include position change time if there are more than 1 positions
    lag = z*2*80 + 4*200 + 840
else
    lag = z*2*80 + 4*200
end
lagperwaitfr = lag/(z*2)
waitperfr = 30000/z-lagperwaitfr*1.2

lasertime = 30000+E*z
uMantime = lasertime*2*P





% %% Multi pos, z scan
% z = 5;
% clc
% P = 12;
% E = 11000;
% if P>1 % include position change time if there are more than 1 positions
%     lag = z*2*80 + 4*200 + 840
% else
%     lag = z*2*80 + 4*200
% end
% lagperwaitfr = lag/(z*2)
% waitperfr = 30000/z-lagperwaitfr*1.2
% 
% lasertime = 30000+E*z
% uMantime = lasertime*2*P
% 
% 
% %% dual timelapse, 1 pos
% z = 1;
% P = 2;
% E = 5000;
% if P>1 % include position change time if there are more than 1 positions
%     lag = z*2*80 + 4*200 + 840
% else
%     lag = z*2*80 + 4*200
% end
% lagperwaitfr = lag/(z*2)
% waitperfr = 30000/z-lagperwaitfr*1.2
% 
% lasertime = 30000+E*z
% uMantime = lasertime*2*P