function[X_sel]=PLSvar_sel(X,Y1,Y2,MaxFac,Type,btnr,btmd,VIPt);

% Variable selection based on bootstrapped-VIP scores calculated from 
% different PLS models 
%
%
% INPUT
% X        Array of independant variables, which is a two way of three way 
%          matrix. Mode 1, Subject; Mode 2, Metabolite; Mode 3, Time.
% Y1       Array of dependent variables representing group infomation.
%          E.g, Samples from intervention and control group are labelled 
%               with 1 and -1 respectively.
% Y2       Array of dependant variables representing time response
%          information.
%          E.g, Samples from response and non-response class are labelled
%               with 10 and 1 respectively.
% MaxFac   Maximal number of components
% Type     Type of model to use:
%          1 bilinear X and group dummy Y
%          2 bilinear X and time-response dummy Y
%          3 bilinear X and group×time-response dummy Y
%          4 trilinear X and group dummy Y
%          5 trilinear X and group/time-response dummy Y
% btnr     Number of bootstrap datasets
% btmd     Resampling method for bootstrap. 1. Balanced resampling 2.
%          Balanced resampling within individual groups. 
% VIPt     VIP threshold
%
% OUTPUT
% X_sel.data    Array of selected variables
% X_sel.index   Index of selected variables
% X_sel.vip     Mean VIP for each variable from bootstrapping

Xorig=X;
switch Type
    case 1
        X=autoscaling(X,2);
        Y=Y1;
    case 2
        X=autoscaling(X,2);
        Y=Y2;
    case 3
        X=autoscaling(X,2);
        Y=Y1.*Y2;
    case 4
        for i=1:size(X,2)
            for j=1:size(X,3)
                X(:,i,j)=autoscaling(X(:,i,j),2);
            end
        end
        Y=Y1;
    case 5
        for i=1:size(X,2)
            for j=1:size(X,3)
                X(:,i,j)=autoscaling(X(:,i,j),2);
            end
        end
        Y=Y1*Y2;
        for i=1:size(Y,1)
            if Y(i,1)==-1
                Y(i,:)=repmat(-1,[1,size(Y,2)]);
            end
        end
end

% Cross validation to choose number of component

[XValResult,Model]=CV(X,Y,MaxFac,1,1);
rms=XValResult.RMSEP;
difference=diff(rms);
ind1=find(difference>0);
ind2=find(abs(difference)<0.01);
if ~isempty([ind1 ind2])
    Fac=min([ind1 ind2]);
else Fac=MaxFac;
end

% Calculate bootrapped VIP
[Vsel]=BTVip(X,Y,Y1,Fac,btnr,btmd,VIPt);
vip=Vsel.mvip;
sel=Vsel.sel;

% Variable selection
[B I]=sort(vip,'descend');
sel=sel(I);
X=Xorig(:,I);
vip=vip(:,I);
index=[1:1:size(X,2)];
index=index(:,I);
if ismember(Type,[1 2 3])
   X_sel.data=Xorig(:,sel);
else
   X_sel.data=Xorig(:,:,sel);
end
X_sel.index=index(1,sel);
X_sel.vip=vip(1,sel);
end


function [Vsel]=BTVip(Xc,Yc,Y1,Fac,btnr,btmd,VIPt);

Ic=size(Xc,1);
DimXc=size(Xc);
% Resampling
switch btmd
    case 1 % choose random samples for bootstrap
        n=[1:1:Ic]';
        L=repmat(n,1,btnr);
        index=[Ic*btnr:-1:2];
        for k=1:(Ic*btnr-1)
            i=index(k);
            j=randi([1,i]);
            L([i j])=L([j i]);
        end
    case 2 % choose random samples from different classes for bootstrap
        ugroup=unique(Y1);
        n=[1:1:Ic]';
        Ind1=n(find(Y1==ugroup(1)));
        Ind2=n(find(Y1==ugroup(2)));
        I1=size(Ind1,1);
        I2=size(Ind2,1);
        
        n1=[1:1:I1]';
        L1=repmat(n1,1,btnr);
        index=[I1*btnr:-1:2];
        for k=1:(I1*btnr-1)
            i=index(k);
            j=randi([1,i]);
            L1([i j])=L1([j i]);
        end
        for i=1:I1
            L1(find(L1==i))=Ind1(i);
        end
        
        n2=[1:1:I2]';
        L2=repmat(n2,1,btnr);
        index=[I2*btnr:-1:2];
        for k=1:(I2*btnr-1)
            i=index(k);
            j=randi([1,i]);
            L2([i j])=L2([j i]);
        end
        for i=1:I2
            L2(find(L2==i))=Ind2(i);
        end
        L=[L1;L2];
end

% Build model and calculate VIP
for i=1:btnr
    Xr=Xc(L(:,i),:);
    Yr=Yc(L(:,i),:);
    Xr=reshape(Xr,DimXc);    
    [Xfactors,Yfactors,~,B,~,~,~,~] = npls(Xr,Yr,Fac,NaN);
    w=Xfactors{1,end}(:,1);
    wabs(i,:)=abs(w);
    model.Xfactors=Xfactors;
    model.Yfactors=Yfactors;
    model.B=B;
    VIP=VIPnway(model, 'nway', 0);
    vip(i,:)=VIP{1,end};
end

% Varibale selection based on VIP
mvip=mean(vip);
svip=std(vip);
lower=mvip-svip;

sel=false(1,DimXc(end));
for i=1:size(vip,2)
    if lower(i)>=VIPt
        sel(i)=true;
    end
end
Vsel.mvip=mvip;
Vsel.svip=svip;
Vsel.lower=lower;
Vsel.sel=sel;

% Variable selection based on w
% mwabs=mean(wabs);
% swabs=std(wabs);
% lowerw=mwabs-swabs;
% 
% selw=false(1,DimXc(end));
% for i=1:size(vip,2)
%     if lowerw(i)>=0
%         selw(i)=true;
%     end
% end
% Vsel.mwabs=mwabs;
% Vsel.swabs=swabs;
% Vsel.lowerw=lowerw;
% Vsel.selw=selw;
end

function [x1sc] = autoscaling(x1,n)
% n=1 mean center
% n= others autoscaling
if n==1
    % Declaration of variables and constants
    [d1rows]=size(x1,1);
    
    % Calculation of mean from X
    xmean=mean(x1);
    
    % Calculation of centered X data
    
    x1sc=x1-(ones(1,d1rows)'*xmean);
    
else
    % Declaration of variables and constants
    [d1rows]=size(x1,1);
    
    % Calculation of std from X
    xstd=std(x1);
    
    % Calculation of mean from X
    xmean=mean(x1);
    
    % Calculation of autoscaled X data
    zeros_x1std=find((ones(d1rows,1)*xstd)==0);
    
    x1sc=(x1-(ones(1,d1rows)'*xmean))./(ones(1,d1rows)'*xstd);
    
    x1sc(zeros_x1std)=0;
end
end






