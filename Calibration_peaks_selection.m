%Script that needs a current Matlab figure opened with the plotted
%calibration spectrum in psuedofrequency [mum^{-1}].
%This script asks the user to type the number of peaks that will be used
%for the calibration. Then the script allows the user to manually select
%the peaks (if the selection of a point goes wrong the user can answer 0 to
%the question if the point has to be kept and reselect the point as many
%times as he wants). For each selected peak the script asks the user to
%type the correspondent calibration value (e.g. in frequency [THz],
%wavelength [nm], wavenumber [cm-1] ect...). At the end the scrip performs a
%Gaussian fit in the around of the points selected by the user: the fit is
%performed with a precision that the user can choose at the beginning.
%According to the fit the exact pseudofrequency is determined. The script
%plots the spectrum showing the performed fits around the selected peaks and
%returns a vector PSEUDOFR with the exact values of pseudofrequency peaks
%and a vector CALIBRATION with the correspondent calibration values.

figure(gcf); %consider the current open figure for the spectrum
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes

pseudofrequency=(get(dataObjs,'XData'));
spectrum=(get(dataObjs,'YData'));

%oppure
% pseudofrequency=cell2mat(get(dataObjs,'XData'));
% spectrum=cell2mat(get(dataObjs,'YData'));

%oppure
% DD=dataObjs{2};
% pseudofrequency=cell2mat(get(DD,'XData'));
% spectrum=cell2mat(get(DD,'YData'));

localMaxima=islocalmax(spectrum); %find all the spectrum peaks

if sum(localMaxima)==0 %if there are not peaks abort
    fprintf('\nThere are not peaks in this spectrum!\n')
    return;
end

prompt={'Select the number of peaks for calibration'};
name='Select number of peaks';
numlines=1;
defaultanswer={num2str(1)};

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

answer=inputdlg(prompt,name,numlines,defaultanswer,options);

num_peaks=str2double(cell2mat(answer));
pos=zeros(num_peaks,1);

PSEUDOFR=zeros(num_peaks,1); %y-axis in calibration (pseudofrequency)
CALIBRATION=zeros(num_peaks,1); %x-axis in calibration (real frequency, wavelength or wavenumber)

precision=pseudofrequency(2)-pseudofrequency(1);

prompt={['Current pseudofrequency precision: ' sprintf('%0.1e',precision) ' \mum^{-1}). Select pseudofrequency calibration precision x (1e-x \mum^{-1})']};
name='Pseudofrequency precision';
numlines=1;
defaultanswer={num2str(7)};

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

answer=inputdlg(prompt,name,numlines,defaultanswer,options);

precision=str2double(cell2mat(answer));

uiwait(msgbox(['Select ' sprintf('%d',num_peaks) ' points. Before each selection ENTER to show the cross.'],'warn'));

index=1;
while index<=num_peaks
    pause;
    [x,~]=ginput(1); %user selected point

    %procedure to find the point in the pseudofrequencies axis closer to x
    array=pseudofrequency;
    array=abs(array-x);
    [~,pos(index)]=min(array); %val: position of x (selected by the user) in psueudofrequencies array

    pos_r=pos(index);
    pos_l=pos(index);

    %find the local maximum closest to pos(index) by exploring on the left (pos_l) and on
    %the right (pos_r)
    while pos_r<=length(localMaxima) || pos_l>0
        if localMaxima(pos_r)==1
            pos(index)=pos_r;
            break;
        end

        if localMaxima(pos_l)==1
            pos(index)=pos_l;
            break;
        end

        if pos_r<length(localMaxima)
            pos_r=pos_r+1;
        end 

        if pos_l>1
            pos_l=pos_l-1;
        end

    end

    %now pos is the position of the peak
    figure(gcf)
    hold on;
    h=scatter(pseudofrequency(pos(index)),spectrum(pos(index)),'o','filled','k');
    
    prompt={'Do you want this point?'};
    name='Do you want this point? [1:yes, 0:no]';
    numlines=1;
    defaultanswer={num2str(1)};

    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';

    answer=inputdlg(prompt,name,numlines,defaultanswer,options);

    risposta=str2double(cell2mat(answer));
    
    if risposta==1
        %the program ask the user to associate the calibration value
        prompt={'Associated peak value'};
        name='Associate a calibration value to the peak';
        numlines=1;
        defaultanswer={num2str(1)};

        options.Resize='on';
        options.WindowStyle='normal';
        options.Interpreter='tex';

        answer=inputdlg(prompt,name,numlines,defaultanswer,options);

        CALIBRATION(index)=str2double(cell2mat(answer));
        
        
        index=index+1;
    else
        delete(h);
    end

end

temp=diff(spectrum)./diff(pseudofrequency);
derivative=zeros(size(spectrum));
derivative(1:end-1)=temp;
derivative(end)=derivative(end-1);

hh=figure;
plot(pseudofrequency,spectrum,'Linewidth',2,'Color','k');

for index=1:num_peaks
    
    %find the extension of the peak slope on the left
    if pos(index)>=3 %there have to be 3 points (2 points in the derivative)
        jj=0;
        while (derivative(pos(index)-1-jj)>0)==(derivative(pos(index)-2-jj)>0)
            jj=jj+1;
            if pos(index)-2-jj<=0
                break;
            end
        end
    end
    
    %find the extension of the peak slope on the right
    if pos(index)<=length(spectrum)-1 %there have to be 1 points (1 points in the derivative)
        kk=0;
        while (derivative(pos(index)+kk)>0)==(derivative(pos(index)+1+kk)>0)
            kk=kk+1;
            if pos(index)+1+kk>length(spectrum)
                break;
            end
        end
    end
    
    xx=pseudofrequency(pos(index)-1-jj:pos(index)+1+kk);
    yy=spectrum(pos(index)-1-jj:pos(index)+1+kk);
    
    f=fit(xx.',yy.','gauss1'); %gaussian fit
    %extrapolated fit parameters
    a1=f.a1;
    b1=f.b1;
    c1=f.c1;
    
    n_points=(max(xx)-min(xx)).*(10^(precision));
    xxx=linspace(min(xx),max(xx),n_points);
    Gaussian_funct=a1.*exp(-((xxx-b1)./c1).^2); %Gaussian Function
    
    [~,Max_Pos]=max(Gaussian_funct); %find Gaussian maximum
    PSEUDOFR(index)=xxx(Max_Pos); %PSEUDOFREQUENCY OF PEAK
    
    figure(hh);
    hold on;
    plot(xxx,Gaussian_funct,'Linewidth',2,'Color','r');
    hold on;
    scatter(PSEUDOFR(index),Gaussian_funct(Max_Pos),'o','filled','r');
    label={['p.f.=' sprintf('%e',PSEUDOFR(index))], ['value=' sprintf('%d',CALIBRATION(index))]};
    text(PSEUDOFR(index),Gaussian_funct(Max_Pos),label,'VerticalAlignment','bottom','HorizontalAlignment','right')
    
end

clear localMaxima prompt name numlines defaultanswer answer num_peaks 
clear pos index x array pos_r pos_l h risposta temp derivative jj kk 
clear xx yy f a1 b1 c1 n_points xxx Gaussian_funct Max_Pos
clear axesObjs dataObjs hh options precision label
clear pseudofrequency spectrum