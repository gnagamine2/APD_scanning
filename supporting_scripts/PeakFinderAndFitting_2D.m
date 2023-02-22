%Created by Gabriel Nagamine
%Last updated 14-09-22

clear all


%Enter folder containing .ASC files from EMCCD for analisys. Examples of
%files can be found on \\nas12.ethz.ch\mavt_ipe_omel_s1\temp\Gabriel\Data\22-09-08_SingleParticle_AP-4-17_LargeAreaConcSample\FilesForPeakAnalysis\

DataDirectory='\\nas12.ethz.ch\mavt_ipe_omel_s1\temp\Gabriel\Data\22-09-08_SingleParticle_AP-4-17_LargeAreaConcSample\FilesForPeakAnalysis\';
AuxScriptDirectory='C:\Users\gabri\Documents\GitHub\2d-peak-finding-and-fitting\'; %Write path to the folder containing extreme and fit gaussian auxiliary scripts

%Enter parameters for analysis.

MinWavelength=510; %for core/alloy MSNC it was 506-552
MaxWavelength=590; %for core/shell QDs it was 580-625
MinYpix=350; %Insert the Y region of the camera that the data should be analysed
MaxYpix=650;
PeakAmplitude=.1*10^6; %Use figures to tune these parameters
MinPeakDistance=5;

%Initiating variables


RawData = [];
DataAux=[];
WAux=[];
Data=[];
widths=[];
centers=[];
errors=[];
Allwidths=[];
Allcenters=[];
Allerrors=[];
Wintegral=[];
AllSpec=[];


%Reading asc files

cd (DataDirectory);
Files=dir('*.*');
Files(1,:)=[]; %Deleting first lines
Files(1,:)=[];

RawData = importdata(Files(1).name);
WAux=RawData(RawData(:,1)>MinWavelength & RawData(:,1)<MaxWavelength, 1);

for i=1:length(Files)
    
    fprintf('Reading file %i/%i \n', i, length(Files))
    RawData = importdata(Files(i).name);
    DataAux=RawData(RawData(:,1)>MinWavelength & RawData(:,1)<MaxWavelength, MinYpix:MaxYpix);
    DataAux=[WAux DataAux];
    Data(:,:,i)=DataAux;
    
end

%% 

%Change to the folder containing extreme and fit gaussian auxiliary scripts
cd (AuxScriptDirectory);

%Processing data

for i=1:length(Files)
    
    fprintf('Processing file %i/%i \n', i, length(Files))
    
    %Integrating data on the wavelength to find the pixel position of the peaks
    Wintegral(:,2,i)=sum(Data(:,:,i));
    Wintegral(1,:,i)=NaN;
    for k=1:length(Wintegral)
        Wintegral(k,1,i)=k;
    end
    
    %Finding the position of the peaks on Y
    [pks,locs,w,p] =findpeaks(Wintegral(:,2, i),'MinPeakProminence',PeakAmplitude,'MinPeakDistance',MinPeakDistance);
    findpeaks(Wintegral(:,2, i),'MinPeakProminence',PeakAmplitude,'MinPeakDistance',MinPeakDistance);
    x=Data(:,1,i);
    
    %figure;%Shows the found peaks after integrating the camera image across the wavelength. Use this figures to tunePeakAmplitude and MinPeakDistance
    
    for j=1:length(locs)
        
        %Creating single spectrum found on Y (locs) position
        SingEmSpec=sum(Data(:,[locs(j,1)-3:locs(j,1)+3]), 2);
        AllSpec(:,j,i)=SingEmSpec;

        %The fitting function has a form of a*exp(-((x-b)/c)^2)+d*x+e
        EmissionEnergy=1240./Data(:,1,i);

        a=7e+04; b=5e+02; c=-12; d=-64; e=1e+05;
        f=a*exp(-((x-b)./c).^2)+d*x+e;
        fn1=SingEmSpec;
        [F1 X1 err it]=fitgaussian(x,fn1);
        [F2 X2]=fitgaussian(EmissionEnergy,F1);

        widths(j,1,i)=X2(1, 3);
        centers(j,1,i)=X2(1, 2);
        errors(j,1,i)=err;
        
        Allwidths=[Allwidths;widths(j,1,i)];
        Allcenters=[Allcenters;centers(j,1,i)];
        Allerrors=[Allerrors;errors(j,1,i)];

        %plot(EmissionEnergy,fn1,'.',EmissionEnergy,F2,'r');
    end

end




%% Generating emission widths, centers, FWHMs
AfterDeletedAllwidths=[];
AfterDeletedAllcenters=[];
AfterDeletedAllwidths(:,1)=Allwidths;
AfterDeletedAllcenters(:,1)=Allcenters;


for j=length(Allwidths):-1:1
    if (Allwidths(j,1)<0.1|Allwidths(j,1)>0.3)
        AfterDeletedAllwidths(j,:)=[];
        AfterDeletedAllcenters(j,:)=[];
        fprintf('Row deleted %i \n', j);
    end

end

AfterDeletedAllfwhm(:,1) = 2*sqrt(log(2))*AfterDeletedAllwidths;



