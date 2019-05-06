function PLOTvar_sel(X_two_way,Y1,Ytime,nr);

% Plot temporal profiles for selected variables
%
% INPUT
% X_two_way    Array of independant variables, which is a two way matrix. 
% Y1           Array of dependent variables representing group infomation.
%              E.g, Samples from intervention and control group are labelled 
%                   with 1 and -1 respectively.
% Ytime        Array of dependant variables representing time
% nr           The index of selected variable resulted from X_sel.index. Input one number a time.     
%
% OUTPUT
% Figure   Temporal profiles for selected variables


% Setting
set(gca,'color',[1 1 1]);
set(gcf,'color',[1 1 1]);
my_colour2 = [148 150 161] ./ 255;
my_colour1 = [124 129 255] ./ 255;

% Input
data=X_two_way(:,nr);
timepoints=sort(unique(Ytime))';

% Processing
datatemp=data;
data1=datatemp(Y1==1,:);
data2=datatemp(Y1==-1,:);
Ytime1=Ytime(Y1==1,:);
Ytime2=Ytime(Y1==-1,:);
for j=1:size(timepoints,2)
    data1t=data1(Ytime1==timepoints(j),:);
    data2t=data2(Ytime2==timepoints(j),:);
    y1(1,j)=mean(data1t);
    e1(1,j)=std(data1t)/sqrt(length(data1t));
    y2(1,j)=mean(data2t);
    e2(1,j)=std(data2t)/sqrt(length(data2t));
end
    
% Plotting
hold on;
xlim([0 max(timepoints)+1]);
set(gca,'XTick',timepoints);
xlabel('Time')
ylabel('Intensity')
p1=plot(timepoints,y1(1,:),'-o');
h1 = errorbar(timepoints,y1(1,:),e1(1,:));
set(p1,...
    'color',my_colour1,...
    'MarkerFaceColor',my_colour1,...
    'MarkerEdgeColor',my_colour1,...
    'MarkerSize',4,...
    'linewidth',1.5)
set(h1,...
    'color',my_colour1,...
    'linewidth',1)

p2=plot(timepoints,y2(1,:),'-^');
h2 = errorbar(timepoints,y2(1,:),e2(1,:));
set(p2,...
    'color',my_colour2,...
    'MarkerFaceColor',my_colour2,...
    'MarkerEdgeColor',my_colour2,...
    'MarkerSize',4,...
    'linewidth',1.5)
set(h2,...
    'color',my_colour2,...
    'linewidth',1)
end
