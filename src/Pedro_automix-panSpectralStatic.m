% Static Spectral Panning of n related sound sources
% Only works with mono files of the same length
% Decisions are based on the probability mass function of each bin
% Placements are made according to scaled roots of chebysheb polynomials.

clear all;
Song_name           = uigetdir;

% Get all files in a directory
a                   = dir([Song_name,'/*.wav']);         




% Populate an array with files in directory
for i= 1:length(dir([Song_name,'/*.wav']))
        txt{i}      = a(i).name;
end

% guarantee all tracks have the same size (ProTools export sometimes adds
% one sample for stereo tracks)
myLen = intmax;
for i = 1:length(txt)
    [x, Fs]         = audioread([Song_name, '/', txt{i}]);
    if size(x,1) < myLen
        myLen       = size(x,1);
    end
    
end
% guarantee all tracks are monophonic and create a matrix with their raw
% values on each column.
for i = 1:length(txt)
    [x, Fs]         = audioread([Song_name, '/', txt{i}]);
    x               = mean(x,2);
    waveforms(:,i)  = x(1:myLen);
end

% Fourier transform specifics
wLen                = 2^10;
hopSize             = wLen/4;
win                 = hanning(wLen, 'periodic');
halfWin             = wLen/2;
k                   = (1:halfWin+1)';
nrWin               = floor((myLen - wLen) / hopSize);

% Create spectral multi-dimensional matrix   
grain               = zeros(wLen,1);
spectrum            = zeros(wLen,1);
magSpectrum         = zeros(wLen,1);
sqMagSpectrum       = zeros(wLen,1);
spectra             = zeros(nrWin+1, wLen, length(txt));
for i = 1:length(txt)
    fprintf('Analyzing track: %s\n',txt{i}); 
    x                       = waveforms(:,i);
    for j = 0:nrWin
        index               = j * hopSize;
        grain               = x(index+1 : index+wLen).* win;
        spectrum            = fft(grain)/halfWin;
        magSpectrum         = abs(spectrum);
        sqMagSpectrum       = magSpectrum.^2;
        spectra(j+1,:,i)    = magSpectrum;
    end  
end

clear i x hopSize windowLength myWindow halfWindow k;
clear j index grain spectrum magnitudeSpectrum squaredMagnitude;
 
% Probability mass functions for spectral bins
avgSpec = nanmean(spectra); % this is something that cannot happen for dynamic

% Find position vectors for each input at each k bin
for i =1:size(avgSpec,2)
    kVector(1,:)                        = avgSpec(1,i,:);
    kVector(2,:)                        = 1:length(kVector);
    kVector                             = sortrows(kVector',1)';
    kVector(1, 1:2:size(kVector,2))     = -kVector(1, 1:2:size(kVector,2));
    kVector                             = sortrows(kVector',1)';
    n                                   = size(avgSpec,3);
    chebyshev                           = [(pi/(2*n)):(pi/n):pi];
    chebyshev                           = sort( cos(chebyshev) );
    chebyshev                           = chebyshev * i / size(avgSpec,2); % defines available positions (-1 to 1)
    positions                           = (chebyshev + 1).* (pi/4); % rescales positions to 0 to pi/2 so as to use sine and cosine
    kVector(3,:)                        = positions;   
    kVector                             = sortrows(kVector', 2)';
    positionVector(i,:)                 = kVector(3,:);
end
 
clear kVector n chebyshev positions;
 
n1              = 256;                          
n2              = n1;                           
s_win           = 1024;                         
w1              = hanning(s_win, 'periodic');   
w2              = w1;   
hs_win          = s_win/2;                      
coef            = sqrt(2)/2; 
a1              = 2*sin((0:hs_win)/s_win*200);
theta           = min(1, max(-1, a1));      
theta2          = theta.' * pi/4;        
theta3          = [theta2(1:hs_win+1) ; flipud(theta2(1:hs_win-1))];


% Read the mono files and turn them into stereo
for i =1:length(txt)
    fprintf('Reading track: %s\n',txt{i});
    x           = waveforms(:,i);
    
    L           = length(x);                    
    x           = [zeros(s_win,1); x; zeros(s_win - mod(L,n1),1)]; % / max(abs(x));
    x_out       = zeros (length(x),2);                               
    
    tic
    pin = 0;
    pout = 0;
    pend = length(x) - s_win;

    while pin < pend
        grain   = x(pin+1: pin+s_win).* w1;          
        f       = fft(grain);
        % Ver bem porque ? que os tamanhos s?o diferentes
        ftL     = coef * f.* (cos(positionVector(:,i)));
        ftR     = coef * f.* (sin(positionVector(:,i)));
        grainL  = (real(ifft(ftL))).* w2;
        grainR  = (real(ifft(ftR))).* w2;
    
        x_out(pout+1: pout+s_win, 1) = x_out(pout+1: pout+s_win, 1) + grainL;
        x_out(pout+1: pout+s_win, 2) = x_out(pout+1: pout+s_win, 2) + grainR;
    
        pin     = pin + n1;
        pout    = pout + n2;
    end
    toc
    
    fprintf('Writing track: %s\n',txt{i});
    x_out       = x_out(s_win+1:s_win+L,:); % / max(max(abs(x_out)));
    %wavwrite(x_out, Fs,['music/newSP/', txt{i}]);
    waveforms_out{i} = x_out;
end

y = zeros(size(waveforms_out{1}));
for i = 1:length(waveforms_out)
    y = y + waveforms_out{i};
end
audiowrite('perdromix.wav',y,Fs)

plot(positionVector(:,1),'r')       % bass
hold on;

plot(positionVector(:,2),'b')
plot(positionVector(:,3),'g')
% plot(positionVector(:,4),'c')
% plot(positionVector(:,5),'m')
% plot(positionVector(:,6),'k')       % snare
% plot(positionVector(:,7),'k')
%plot(positionVector(:,8),'w')
%plot(positionVector(:,9),'b')     % guitar
%plot(positionVector(:,10),'b--')
%plot(positionVector(:,11),'g--')
%plot(positionVector(:,12),'g')    % shaker
%plot(positionVector(:,13),'m--')
%plot(positionVector(:,14),'y--')
hold off;
% 
