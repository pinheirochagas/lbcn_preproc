function wave = noiseFiltData(globalVar, wave)
% Filtering 60 Hz line noise
% Dependencies:notch.m .eeglab toolbox
% Writen by Mohammad Dastjerdi, Parvizi Lab, Stanford
% Revision date SEP,2009



% filtering 60 Hz
if strcmp(globalVar.center, 'Stanford')
    [wave]= notch(wave, globalVar.iEEG_rate, 59, 61,1);
    [wave]= notch(wave, globalVar.iEEG_rate, 118,122,1); % Second harmonic of 60
    [wave]= notch(wave, globalVar.iEEG_rate, 178,182,1); % Third harmonic of 60
%     
%     d = designfilt('bandstopiir','FilterOrder',2, ...
%         'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
%         'DesignMethod','butter','SampleRate',globalVar.iEEG_rate);
%     wave = filtfilt(d,double(wave));
%     d = designfilt('bandstopiir','FilterOrder',2, ...
%         'HalfPowerFrequency1',118,'HalfPowerFrequency2',122, ...
%         'DesignMethod','butter','SampleRate',globalVar.iEEG_rate);
%     wave = filtfilt(d,double(wave));
%     d = designfilt('bandstopiir','FilterOrder',2, ...
%         'HalfPowerFrequency1',178,'HalfPowerFrequency2',182, ...
%         'DesignMethod','butter','SampleRate',globalVar.iEEG_rate);
%     wave = filtfilt(d,double(wave));
%     

    
elseif strcmp(globalVar.center, 'China')
    [wave]= notch(wave, globalVar.iEEG_rate, 49, 51,1);
    [wave]= notch(wave, globalVar.iEEG_rate, 98,102,1); % Second harmonic of 60
    [wave]= notch(wave, globalVar.iEEG_rate, 148,152,1); % Third harmonic of 60
else
end


end