function  [VIP]=VIPnway(model, tool, plots);

% Calculate Variable Inluence in Projection (VIP) from a NPLS regression
% model.
% It works for X n-way and Y up to two-way and it assume samples are in the first
% mode.
% VIP estimate the importance of each variable in the projection used in a NPLS model and it could be used for
% variable selection. A variable with a VIP value close to or greater than
% 1 (one) can be considered important in given model. 
% However, this is just an empirical rule and may need to be adapted to the
% nature of the analyzed data set.
%
% %I/O: VIP_nway = VIPnway(model, tool,plots)
%
% INPUTS:
% 1. model
%    model = a model structure calcualted by NPLS of PLS_Toolbox
%%%% e.g. model= npls(X,Y,nlv);
%    OR if you are using npls from the nway3.1 toolbox
%%%% e.g. [Xfactors,Yfactors,Core,B] = npls(X,Y,nlv);
%    you should assemble model befor erunning this code in this way:
%      model.Xfactors=Xfactors;
%      model.Yfactors=Yfactors;
%      model.B=B;
% 2. tool
%      tool= flag for telling if model come from PLS_Toolbox  'plst' or
%      nway3.1 toolbox 'nway'
%
% OPTIONAL INPUT:
% 3. plots
%   plots = flag for plots:  0  no plots in output; 1   VIP^2 plots in
%   output
%  
% OUTPUTS: VIPscores a cell array with VIP, for each Y, from mode 2 to mode N in each cell
%    VIPscores.partial[y's, mode]  
%   a cell array of dimension: number of y times number of modes
%   each row cell holds VIP for the corresponding y variables,
%   the columns refer to the different variables/conditions modes of the array.
%
%    VIPscores.all[calculated as a summarized score for all Y, mode2 to n]
%
%%%%  in case of only one y VIPscores.partial=VIPscores.all

%%  $ Version 1 $ Feb 2014 $ Marina Cocchi $ Not compiled $
%%   please cite:  
%%   S. Favilla C. Durante,M. Li Vigni,M. Cocchi, 
%%  Assessing feature relevance in NPLS models by VIP
%%   Chemom. Intell. Lab. Syst. 2013 (129) 76-86.
%%


% if ~exist('plots','var');
%       plots=0;
% end

switch tool
    case {'plst'}
        dim_arX=size(model.loads,1); % order of the array
        
        T=model.loads{1,1};
        nlv=size(T,2);
        
        Q=model.loads{2,2};
        B=model.detail.bin; 
        
        dim_y=size(Q,1); 
        
        for j=2:dim_arX;
            w{j}=model.loads{j,1};  % NPLS weights of mode j: J(or K,..) x LVs
            dim_w(j)=size(w{j},1); 
        end
        
    case {'nway'}
        dim_arX=length(model.Xfactors);
        T=model.Xfactors{1};
        nlv=size(T,2);
        for j=2:dim_arX;
            w{j}=model.Xfactors{j};  % NPLS weights of mode j: J(or K,..) x LVs
            dim_w(j)=size(w{j},1);  
        end
        Q=model.Yfactors{2};
        B=model.B;
        dim_y=size(Q,1);
end
%%-------- maybe not needed
 for j=2:dim_arX;
  wnorm{j}=(w{j}*diag(1./sqrt(sum(w{2}.^2,1))));
 end
%%---------

% start calculation for VIP, first part is used, for each Y, to assess 
% Y variance explained by Factor i with respect to Y variance explained
% by Factor from 1 to nlv

yy=T*B; %[I nLV]x[nLV nLV] = [I nLV]
  
 for i=1:nlv; 
  yyy(:,:,i)=yy(:,i)*Q(:,i)'; %[I 1]x [1 x Y] -> [I Y nLV] 
 end 

syyy=sum(yyy.^2); % sum over I -> [1 Y nLV]
syyy_shdim=shiftdim(syyy); %[Y nLV] variance of each y for each LV

 if dim_y > 1;
  syyy_shdim_perm=permute(syyy_shdim, [2 1]); % invert first second dimension [nLV Y]
 else
  syyy_shdim_perm= syyy_shdim;
 end
  % ssyyy=(sum(syyy_shdim_perm)); % Bug. Changed to the line below 
  ssyyy=(sum(syyy_shdim_perm,1)); % [1 Y] total variance for each y

 for k=1:dim_y;
  var_y(:,k)=syyy_shdim_perm(:,k)./ssyyy(k)'; % [LV 1]./[1] -> [LV Y]
 end

 for j=2:dim_arX;
  for k=1:dim_y;
   VIP_y{k,j-1} = dim_w(j)*(w{:,j}.^2*var_y(:,k));
  end
 end

%%%% VIP for all Y altogether

 for i=1:nlv; syyy2(i)=sum(sum(yyy(:,:,i).^2));end
 
 for j=2:dim_arX;
  VIP{j} = (dim_w(j)*(wnorm{j}.^2*(syyy2./sum(sum(yy.^2)))'));
 end

%  VIPscores.partial=VIP_y;
%  VIPscores.all=VIP;
 

%  switch plots
%      case 1
%          for j=2:dim_arX; % mode 2 to n
%              if tool=='plst'
%                  if  ~isempty(model.detail.axisscale{j,1})
%                      axis_mode=model.detail.axisscale{j,1};
%                  else
%                      axis_mode=[1:dim_w(j)];
%                  end
%              elseif tool=='nway'
%                  axis_mode=[1:dim_w(j)];
%              end
%              
%              
%              min_ax=min(axis_mode);
%              max_ax=max(axis_mode);
%              
%              
%              figure;
%              for k=1:dim_y; % number of y
%                  subplot(ceil(sqrt(dim_y)), ceil(sqrt(dim_y)), k);
%                  hh=plot(axis_mode, VIPscores.partial{k,j-1}, 'b','LineWidth',1.5);
%                  hold on
%                  plotline1=hline(1, 'r');
%                  set(plotline1,'LineWidth',2);
%                  if dim_y == 1;
%                      title(['mode' num2str(j) ' vip y' num2str(k) '=Y : ' num2str(nlv) 'LV'], 'FontWeight','bold','FontSize',12);
%                  else
%                      title(['mode' num2str(j) ' vip y' num2str(k) ': ' num2str(nlv) 'LV'], 'FontWeight','bold','FontSize',12);
%                  end
%                  axis([min_ax max_ax 0 3]);
%              end
%              if dim_y > 1
%                  figure;
%                  plot(axis_mode, VIPscores.all{1,j},'k', 'LineWidth',1.5);
%                  hold on
%                  plotline2=hline(1, 'r');
%                  set(plotline2,'LineWidth',2);
%                  title(['mode' num2str(j) ' vip all Y: ' num2str(nlv) 'LV'], 'FontWeight','bold','FontSize',12);
%                  axis([min_ax max_ax 0 3]);
%              end
%          end
%  end

% save VIPscores VIPscores