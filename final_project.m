%%
%Part1: Using Matlab to select the potential interacted proteins with EV71 3D 
%polymerase based on the raw data.

%step1: calculate the spectrum number for each gene from row data
data1=readtable('3DFlag_EV71.txt');
file=data1(:,[1,4,5]);
file=table2cell(file);
file(end+1,1)={0}; %add one line for calculating the spectrum number of last gene 
count=1;
FlagEV71={};
j=0;
for i=2:length(file)
    if cell2mat(file(i-1,1))==cell2mat(file(i,1))
        count=count+1; %spectrum number of one gene
    else
        j=j+1; %one j for one gene
        FlagEV71(j,1:2)=file(i-1,2:3);
        FlagEV71(j,3)=num2cell(count);
        count=1;
    end  
end


data2=readtable('vector_EV71.txt');
file=data2(:,[1,4,5]);
file=table2cell(file);
file(end+1,1)={0}; %add one line for calculating the spectrum number of last gene 
count=1;
vectorEV71={};
j=0;
for i=2:length(file)
    if cell2mat(file(i-1,1))==cell2mat(file(i,1))
        count=count+1; %spectrum number of one gene
    else
        j=j+1; %one j for one gene
        vectorEV71(j,1:2)=file(i-1,2:3);
        vectorEV71(j,3)=num2cell(count);
        count=1;
    end  
end


data3=readtable('3DFlag_noEV71.txt');
file=data3(:,[1,4,5]);
file=table2cell(file);
file(end+1,1)={0}; %add one line for calculating the spectrum number of last gene 
count=1;
FlagnoEV71={};
j=0;
for i=2:length(file)
    if cell2mat(file(i-1,1))==cell2mat(file(i,1))
        count=count+1; %spectrum number of one gene
    else
        j=j+1; %one j for one gene
        FlagnoEV71(j,1:2)=file(i-1,2:3);
        FlagnoEV71(j,3)=num2cell(count);
        count=1;
    end  
end


data4=readtable('vector_noEV71.txt');
file=data4(:,[1,4,5]);
file=table2cell(file);
file(end+1,1)={0}; %add one line for calculating the spectrum number of last gene 
count=1;
vectornoEV71={};
j=0;
for i=2:length(file)
    if cell2mat(file(i-1,1))==cell2mat(file(i,1))
        count=count+1; %spectrum number of one gene
    else
        j=j+1; %one j for one gene
        vectornoEV71(j,1:2)=file(i-1,2:3);
        vectornoEV71(j,3)=num2cell(count);
        count=1;
    end  
end

%%
%step2: For each gene in the FlagEV71 and FlagnoEV71 groups, find the
%corresponding gene in the control groups (vectorEV71 and vectornoEV71,
%respectively).The spectrum numbers of genes in control groups are considered 
%as background.
FlagEV71(:,4)={0};
for i=1:length(FlagEV71)
    geneaces=cell2mat(FlagEV71(i,1));
    for j=1:length(vectorEV71)
        if contains(cell2mat(vectorEV71(j,1)),geneaces) 
            FlagEV71(i,4)=vectorEV71(j,3);
        end
    end
end


FlagnoEV71(:,4)={0};
for i=1:length(FlagnoEV71)
    geneaces=cell2mat(FlagnoEV71(i,1));
    for j=1:length(vectornoEV71)
        if contains(cell2mat(vectornoEV71(j,1)),geneaces) 
            FlagnoEV71(i,4)=vectornoEV71(j,3);
        end
    end
end
    
%%
%step3: Selecting proteins in the samples of 3DFlag_EV71 and 3DFlag_noEV71
%with spectra number>=4, which indicates that those proteins are indeed
%detected by mass spectrum (not random aligned). Next, selecting proteins with
%fold change>=2 comparing with control samples, and those are the final potential 
%interacted proteins with 3D polymerase.

%spectra number>=4
excul=[];
for i=1:length(FlagEV71)
    if cell2mat(FlagEV71(i,3))<4
        excul(end+1)=i;
    end
end
FlagEV71(excul,:)=[];

excul=[];
for i=1:length(FlagnoEV71)
    if cell2mat(FlagnoEV71(i,3))<4
        excul(end+1)=i;
    end
end
FlagnoEV71(excul,:)=[];

%%
%fold change>=2
for i=1:length(FlagEV71)
    if cell2mat(FlagEV71(i,4))==0 %replace 0 by 1 so that they could serve as denominator
        FlagEV71(i,4)={1};
    end
end

for i=1:length(FlagnoEV71)
    if cell2mat(FlagnoEV71(i,4))==0 %replace 0 by 1 so that they could serve as denominator
        FlagnoEV71(i,4)={1};
    end
end
        
potint_FlagEV71={};
j=0; %line number of potint_FlagEV71
for i=1:length(FlagEV71)
    fodchg=cell2mat(FlagEV71(i,3))/cell2mat(FlagEV71(i,4));
    if fodchg>=2
        j=j+1;
        potint_FlagEV71(j,:)=FlagEV71(i,:);
    end
end
writecell(potint_FlagEV71,'potint_FlagEV71.txt');

potint_FlagnoEV71={};
j=0;
for i=1:length(FlagnoEV71)
    fodchg=cell2mat(FlagnoEV71(i,3))/cell2mat(FlagnoEV71(i,4));
    if fodchg>=2
        j=j+1;
        potint_FlagnoEV71(j,:)=FlagnoEV71(i,:);
    end
end
writecell(potint_FlagnoEV71,'potint_FlagnoEV71.txt');

%%
%step4: visualize the data
subplot(1,2,1);
plot(cell2mat(FlagEV71(:,4)),cell2mat(FlagEV71(:,3)),'.','MarkerEdgeColor',...
[0 0.4470 0.7410],'MarkerSize',10);
hold on
plot(cell2mat(potint_FlagEV71(:,4)),cell2mat(potint_FlagEV71(:,3)),'.','MarkerEdgeColor',...
[0.4660 0.6740 0.1880],'MarkerSize',10);
plot(cell2mat(potint_FlagEV71(7,4)),cell2mat(potint_FlagEV71(7,3)),'.','MarkerEdgeColor',...
[0.8500 0.3250 0.0980],'MarkerSize',10);
hold on
plot([0 90],[0 180],'k-');
xlabel('Spectra Number for vector-EV71');
ylabel('Spectra Number for 3DFlag-EV71');
title('3DFlag-EV71');
legend('Background proteins','Potential interacted proteins','Bait proteins');

subplot(1,2,2);
plot(cell2mat(FlagnoEV71(:,4)),cell2mat(FlagnoEV71(:,3)),'.','MarkerEdgeColor',...
[0 0.4470 0.7410],'MarkerSize',10);
hold on
plot(cell2mat(potint_FlagnoEV71(:,4)),cell2mat(potint_FlagnoEV71(:,3)),'.','MarkerEdgeColor',...
[0.4660 0.6740 0.1880],'MarkerSize',10);
plot(cell2mat(potint_FlagnoEV71(4,4)),cell2mat(potint_FlagnoEV71(4,3)),'.','MarkerEdgeColor',...
[0.8500 0.3250 0.0980],'MarkerSize',10);
plot([0 125],[0 250],'k-');
xlabel('Spectra Number for vector-noEV71');
ylabel('Spectra Number for 3DFlag-noEV71');
title('3DFlag-noEV71');
legend('Background proteins','Potential interacted proteins','Bait proteins');

%%
%Another way to visualize the data (fold change)
fodchg1=cell2mat(FlagEV71(:,3))./cell2mat(FlagEV71(:,4));
[~,I]=max(fodchg1);
fodchg1(I)=[];

fodchg2=cell2mat(FlagnoEV71(:,3))./cell2mat(FlagnoEV71(:,4));
[~,I]=max(fodchg2);
fodchg2(I)=[];

subplot(1,2,1);
hold on;
m=length(fodchg1);
for i=1:m
    if fodchg1(i)>=2
       plot(i,fodchg1(i),'r.','MarkerSize',15,'MarkerEdgeColor',...
       [0.8500 0.3250 0.0980]);
    else
       plot(i,fodchg1(i),'b.','MarkerSize',8,'MarkerEdgeColor',...
       [0.3010, 0.7450, 0.9330]);
    end
end
xlim([0,length(FlagEV71)]);
ylabel('flod change');
title('3DFlag-EV71')
plot([0 m],[2 2],'k-','LineWidth',1);
hold off


subplot(1,2,2);
hold on;
m=length(fodchg2);
for i=1:m
    if fodchg2(i)>=2
       plot(i,fodchg2(i),'r.','MarkerSize',15,'MarkerEdgeColor',...
       [0.8500 0.3250 0.0980]);
    else
       plot(i,fodchg2(i),'b.','MarkerSize',8,'MarkerEdgeColor',...
       [0.3010, 0.7450, 0.9330]);
    end
end
xlim([0,length(FlagnoEV71)]);
ylabel('flod change');
title('3DFlag-noEV71')
plot([0 m],[2 2],'k-','LineWidth',1);
hold off


%%
%Part2: Comparing the differential proteins in the two groups (3DFlag_EV71 
% & 3DFlag_noEV71) and display the result by Venn diagram. (i.e. which proteins 
%are included in both groups, which are not).
count=0;
for i=1:length(potint_FlagEV71)
    for j=1:length(potint_FlagnoEV71)
        if contains(cell2mat(potint_FlagEV71(i,1)),cell2mat(potint_FlagnoEV71(j,1)))
            count=count+1;
        end
    end
end

EV71only=length(potint_FlagEV71)-count;
noEV71only=length(potint_FlagnoEV71)-count;

F = struct('Display', 'iter'); %venn function downloaded from internet
[H,S] = venn([EV71only, noEV71only],count,F,'ErrMinMode','ChowRodgers','FaceAlpha', 0.6); 
text(S.ZoneCentroid(1,1)-1, S.ZoneCentroid(1,2), [num2str(EV71only)]);%label venn
text(S.ZoneCentroid(2,1)+0.5, S.ZoneCentroid(2,2), [num2str(noEV71only)]);
text(S.ZoneCentroid(3,1), S.ZoneCentroid(3,2), [num2str(count)]);
text(-1,3.5,['3DFlag-EV71']);
text(3,3.5,['3DFlag-noEV71']);
title('Differential proteins comparation');

%%
%Part3: Molecular function classification.
%Firstly, using the PANTHER website to classify protein functions. Then, analysis
%results in Matlab.
data1=readtable('3DFLAGEV71 PANTHER Molecular Function.csv');
data2=readtable('3DFLAGnoEV71 PANTHER Molecular Function.csv');

subplot(2,1,1);
file=table2cell(data1(:,3));
geneNumber=cell2mat(file);
pie(geneNumber);
labels=table2cell(data1(:,2));
legend(labels);
title('3DFlag-EV71');

subplot(2,1,2);
file=table2cell(data2(:,3));
geneNumber=cell2mat(file);
pie(geneNumber);
labels=table2cell(data2(:,2));
legend(labels);
title('3DFlag-noEV71');

%%
%Part4: Using KOBAS website to perform KEGG pathway enrichment and analysis
%results in Matlab. Results visualized by bubble map.
data=readtable('KEGG.csv');
file=table2cell(data(1:20,[1 4 5 6])); %select the top 20 pathway to display
richfactor=cell2mat(file(:,2))./cell2mat(file(:,3));%ratio of the enriched number and the number of genes have been annotated in this pathway. 
pvalue=-log10(cell2mat(file(:,4)));

figure; hold on;%plot 4D bubble map

x=richfactor;
y_values=[1:20];
y_labels=file(:,1);
sz=cell2mat(file(:,2));%enriched gene number as circle size in bubble map
c=pvalue;%-log10 pvalues as the color map for circles

bubsizes=unique(sz)';%create size legend for bubble map using 'plot'. It's stupid. But there is no other way to do that?
bubsize=[2.7,3.4,4.1,4.8];%the 'MarkerSize' for plot and 'sz' for scatter use different algorithm..., need to alter the value to get the smae size
legentry=cell(size(bubsizes));
for i = 1:length(bubsizes)
   bubleg(i) = plot(0,0,'k.','markersize',bubsize(i)*8.5);
   set(bubleg(i),'visible','off')
   legentry{i} = num2str(bubsizes(i));
end

scatter(x,y_values,30*sz,c,'filled');
xlim([0 0.065]);
ylim([0 21]);
cb=colorbar; title(cb,'-log10(Pvalue)');%colorbar 
grid on; 
colormap cool;
set(gca,'Ytick',y_values,'YtickLabel',y_labels);
xlabel('Rich Factor');
ylabel('Patheay Name');
title('Pathway Enrichment');
lg=legend(legentry); title(lg,'Number'); 