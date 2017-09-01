%% TransitionMatrixCalculator

%% Inputs

%Load in BTC1 and BTC2 data as well as the actual transition Matrix, called
%M. 
%Pe = 1000
load('count1_1000.mat')
load('count2_1000.mat')
load('tau1_1000.mat')
load('tau2_1000.mat')

%Inputs
dec = 3; %Preprocessing data: the decimal place to which to round your arrival times. In this example 3 works best
bn = 10; %number of transition classes. This will estimate a transition matirx of bn^2

%% 
unv1= tau1 ;  
unv2 =tau2 ;
btc1 = unv1;  
btc2 = unv2;


%% Finding PDFs of arrival times 

%PDF1
CDF1 = zeros(length(tau1),1); %Create Cumulative Distribution Function for BTC at x=L
PDF1 = CDF1; %Create Probability Density Function for BTC at x=L
for i = 2:length(tau1);
    a = i;
    CDF1(a) = trapz(tau1(1:a),count1(1:a));
    PDF1(a) = CDF1(a) - CDF1(a-1);
end
PDF1(1) = 1 - sum(PDF1);

%PDF2
CDF2 = zeros(length(tau2),1); %Create Cumulative Distribution Function at x=2L
PDF2 = CDF2; %Create Probability Density Function at x=2L
for i = 2:length(tau2);
    a = i;
    CDF2(a) = trapz(tau2(1:a),count2(1:a));
    PDF2(a) = CDF2(a) - CDF2(a-1);
end
PDF2(1) = 1 - sum(PDF2);

%% Classifying Arrival Times
%Determine the time that separate velocity classes. 
%The ith cutoff is the slowest time in the ith class
cut =zeros(10,1);
for i = 1:9
    a =abs(ones(length(tau1),1)*.1*i-CDF1); 
    c = find(a == min(a));
    cut(i) = tau1(c);
end
cut(10) = max(tau1);

% Round times to the specified dec. This smooths out noise in the measured BTCs 
tau1r= round(tau1,dec); %round time 1 to the nearest dec
tau2r = round(tau2,dec); %round time 2 to the nearest dec

tau1 = unique(tau1r); %find the unique times for BTC1 after rounding
c1r =[]; % rounded count
tau2 = unique(tau2r); %find unique times BTC2 after rounding
c2r = []; %rounded count
pdf1 = []; %rounded pdf
pdf2 = [];
for i = 1:length(tau1)
   a= find(tau1r == tau1(i)); 
   c1r = [c1r; sum(count1(a))];
   pdf1 = [pdf1; sum(PDF1(a))];
end
tau1 = round(tau1,dec);
count1=c1r;
PDF1 =pdf1;

for i = 1:length(tau2)
   a= find(tau2r == tau2(i)); 
   c2r = [c2r; sum(count2(a))];
   pdf2 = [pdf2; sum(PDF2(a))];
end
count2=c2r;
tau2 = round(tau2,dec);
PDF2 = pdf2;



%% Organize Data
%Determine the amount of particles in each bin
bin_size = 1/bn; 
time_mod = tau1; %used to modify the time vector; Note that some times can be 
%contained in two classes (this may be necessary to ensure velocity classes
%have equal probability). time_mod accounts for these times


%Create a new pdf and time vector that is divides times that occur in
%multiple bins

% When a time occurs in more than one bin, we treat the time as separte to
% each bin.
for i=1:bn-1
pdf_tick = 0; %finds where PDFs are divided into next class
for j = 1:length(PDF1)
    pdf_tick = pdf_tick+PDF1(j);
    if pdf_tick > i*bin_size
           PDF1 = [PDF1(1:j-1); PDF1(j)-(pdf_tick-i*bin_size);pdf_tick-i*bin_size; PDF1(j+1:end)];
           time_mod = [time_mod(1:j-1); time_mod(j);time_mod(j);time_mod(j+1:end)];
           break
    elseif pdf_tick == i*bin_size
        break
    end
end
end

a= find(abs(PDF1)<10^-14); %Correct for Matlab rounding errors
time_mod(a) = [];
PDF1(a)=[];

%% Finding Conditional Probabilites
%Make CDF
CDF = PDF1; 
for i = 2:length(PDF1)
    CDF(i) = CDF(i-1) + PDF1(i);
end

%Determine where the bin where each bin ends
bin_cut_off = zeros(bn,1);
for i = 1: bn-1
    bin_cut_off(i) = find((abs(i*bin_size - CDF))<10^-12);
end
bin_cut_off(bn) = length(CDF);

%Find the distribution for each bin
pdf_cut_off = PDF1;
for i = 2:bn
    pdf_cut_off(bin_cut_off(i-1)+1:bin_cut_off(i)) = ...
        pdf_cut_off(bin_cut_off(i-1)+1:bin_cut_off(i))/sum(pdf_cut_off(bin_cut_off(i-1)+1:bin_cut_off(i)));
end
pdf_cut_off(1:bin_cut_off(1))=pdf_cut_off(1:bin_cut_off(1))/sum(pdf_cut_off(1:bin_cut_off(1)));

%Create a vector where each element correpsonds to the bin that each time
%is in
bin_number = zeros(length(time_mod),1); %the bin that corresponds to each time in time_mod
bin_number(1:bin_cut_off(1))=1;
for i = 2:bn
   bin_number(bin_cut_off(i-1)+1:bin_cut_off(i))=i; 
end

%% Construction of Systems of Equations

%%Create a matrix to solve for transitions probs
IM = zeros(length(PDF2),bn*bn); %Probability Contribution Matrix
time_mod = round(time_mod,dec);
for i = 1:length(PDF2) %for each discrete time in PDF2
    for j = 1:length(time_mod) %for each time from PDF1
        a = round(tau2(i) - time_mod(j),dec); %find the difference between all 
        %combos of times in BTC2 and BTC1. This gives a time in the
        %trvel time distribution of cell 2        
        a = find(time_mod == a); % Note cell2 travel times distribution is equal to 
        %the arrival time distribition for BTC1. So we can assign a
        %velocity class for the cell 2 travel times
        
        b = bin_number(j); %the bin we are in
        c = bin_number(a); %the bin we are transitioning to
        %IM is probabilty contribution for each combo
        if length(a) == 1
            IM(i,bn*(b-1)+c) = IM(i,bn*(b-1)+c)+ PDF1(j)*pdf_cut_off(a);
        elseif length(a)==2
            IM(i,bn*(b-1)+c(1)) = IM(i,bn*(b-1)+c(1))+ PDF1(j)*pdf_cut_off(a(1));
            IM(i,bn*(b-1)+c(2)) = IM(i,bn*(b-1)+c(2))+ PDF1(j)*pdf_cut_off(a(2));
        end
    end
end

