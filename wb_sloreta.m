function [JJ2] = wb_sloreta(data,leadfield,alpha)
% EEG source imaging using sLORETA method (based on Dale et al.standardization).
%  case 1: EEG with unknown current density vector.
%  case 2: EEG with known current density vector orientation, unknown
%          amplitude.This case usually corresponds to the inverse problem 
%          when the cortical surface is completely known. Voxels are now 
%          distributed along the cortical surface, and the dipoles at each 
%          voxel have known orientation (perpendicular to the cortical surface).
%
% Input:
%      data: observation EEG signals.Demension is channels X time points,
%            e.g. 60 channels X 1000 time points
%      leadfield: the lead field. The leadfield can be a matrix
%            (channels X sources/dipoles, e.g. 60 channels X 6144 sources/dipoles)
%            or a matlab structure containing leadfield of x,y,z-orientations
%            (e.g. leadfield.X with demension channels X dipoles (x-orientation); 
%            leadfield.Y with demension channels X dipoles (y-orientation); 
%            leadfield.Z with demension channels X dipoles (z-orientation))
%            which is calculated by using the forward theory, based on the
%            electrode montage, head model and equivalent source model. It
%            can also be the output of ft_prepare_leadfield.m (e.g. lf.leadfield,
%            dipoles contain x,y,z-orientations, 60 channels X 6144*3 dipoles)
%            based on real head modal (FEM modal) using FieldTrip.
%      alpha: regularization parameter.  alpha >= 0: use user-defined value;
%             alpha = 'mean': use averaged alpha with tikh regularization (default);
%             alpha = 'timevarying': use time-varying alpha with tikh regularization.
%             The time cost of 'time-varying' alpha is extremely high for EEG
%             time courses. It is used for the situation of less topographies.
% Output:
%      JJ2: The primary (impressed) current density J ¡ÊR£¨Nd X Nt)
%           The demension is  No. of dipoles X No.of time points
%
% Reference:
% R.D. Pascual-Marqui. Standardized low resolution brain electromagnetic tomography
%            (sLORETA): technical details. Methods & Findings in Experimental
%             & Clinical Pharmacology 2002, 24D:5-12.
% Dale AM, Liu AK, Fischl BR, Buckner RL, Belliveau JW, Lewine JD, Halgren E.
%             Dynamic statistical parametric mapping: combining fMRI and 
%             MEG for highresolution imaging of cortical activity. 
%             Neuron 2000, 26: 55-67.
% -------------------------------------------------------------------------
% Written by Li Dong (UESTC, Lidong@uestc.edu.cn)
% $ 2020.12.9
% -------------------------------------------------------------------------

if nargin < 2
  error ('2 inputs are reqiured at least!!!!!');
elseif nargin == 2
  alpha = 'mean';
end

% check alpha
if isnumeric(alpha)
  if isempty(alpha) || alpha < 0
    disp('alpha value is invalid, use averaged alpha with tikh regularization...')
    alpha = 'mean';
  end
end
% check data
if all(~isfinite(data(:))) || isempty(data)
  error('The data contains NaN or Inf, or it is empty!!!');
end

% check and get the leadfield matrix
flag1 = 0; % flag of whether containing 3 xyz oritations in leadfield
if isnumeric(leadfield)
  G = leadfield; % channels X dipoles (a oritation, e.g. normal vector)
elseif isfield(leadfield,'X') && isfield(leadfield,'Y') && isfield(leadfield,'Z')
  G = [leadfield.X,leadfield.Y,leadfield.Z];      
  flag1 = 1; % contains xyz oritations
      % the leadfield matrix (chans X dipoles*3 ), which
      % contains the potential or field distributions on all
      % sensors for the x,y,z-orientations of the dipole.
elseif isfield(leadfield,'lf') % ONLY for FieldTrip
  if isstruct(leadfield.lf)
    try
      Npos = size(leadfield.lf.pos,1);
      m = 1;
      for i = 1:Npos
        if ~isempty(leadfield.lf.leadfield{1,i})
          lf_X(:,m) = leadfield.lf.leadfield{1,i}(:,1); % X orientation of the dipole.
          lf_Y(:,m) = leadfield.lf.leadfield{1,i}(:,2); % Y orientation of the dipole.
          lf_Z(:,m) = leadfield.lf.leadfield{1,i}(:,3); % Z orientation of the dipole.
          m = m + 1;
        end
      end
      G = [lf_X,lf_Y,lf_Z];
      flag1 = 1; % contains xyz oritations
      % the leadfield matrix (chans X dipoles*3 ), which
      % contains the potential or field distributions on all
      % sensors for the x,y,z-orientations of the dipole.
    catch
      error('leadfiled is not calculated by ''ft_prepare_leadfield.m'' (dipoles contain x,y,z-orientations)?');
    end
  end
elseif isfield(leadfield,'leadfield') % ONLY for FieldTrip
  if iscell(leadfield.leadfield)
    try
      Npos = length(leadfield.leadfield);
      m = 1;
      for i = 1:Npos
        if ~isempty(leadfield.leadfield{1,i})
          lf_X(:,m) = leadfield.leadfield{1,i}(:,1); % X orientation of the dipole.
          lf_Y(:,m) = leadfield.leadfield{1,i}(:,2); % Y orientation of the dipole.
          lf_Z(:,m) = leadfield.leadfield{1,i}(:,3); % Z orientation of the dipole.
          m = m + 1;
        end
      end
      G = [lf_X,lf_Y,lf_Z];
      flag1 = 1; % contains xyz oritations
      % the leadfield matrix (chans X dipoles*3), which
      % contains the potential or field distributions on all
      % sensors for the x,y,z-orientations of the dipole.
    catch
      error('leadfiled is not calculated by ''ft_prepare_leadfield.m'' (dipoles contain x,y,z-orientations)?');
    end
  end
end

% check dimension
DIM1 = size(data);
DIM2 = size(G);

if DIM1(1)~=DIM2(1)
  error('The No. of channels in the data is not equal to the demension of channels in the leadfield matrix, please check your data or leadfield!!!');
end
Nc = DIM1(1); % No. of channels
Nt = DIM1(2); % No. of time points OR No. of topographies
Nd = DIM2(2); % No. of dipoles or dipoles*3 correponding to x,y,and z oritations.

% -------------------------------------------------------------------------
% inverse
JJ1 = zeros(Nd,Nt);  % inverse results
H = eye(Nc)-(ones(Nc,1)*ones(Nc,1)')/(ones(Nc,1)'*ones(Nc,1)); % Eq. 7: H = I-11^T/1^T1
K = H * G;           % averaged transforms
phi = H * data;      % average data (re-referenced to AVG)

% regularization parameter
if ischar(alpha)
  if isequal(alpha,'mean') % sloreta with averaged alpha with tikh regularization
    [U,s,~] = csvd(G); % standard-form regularization
    alpha = zeros(1,Nt);
    for k = 1:Nt
      alpha(k) = l_curve(U,s,data(:,k),'tikh');
    end
    alpha = mean(alpha);
    
    % standardiaztion
    T = K'*pinv(K*K' + alpha*H); % Eq. 6 or Eq.11: T = K^TH[HKK^T + alpha*H]+
    S_j = T * K;  % Eq. 17: S_j = K^TH[HKK^T + alpha*H]+
    J1 = S_j * (T * phi);
    
    % Using Dale et al.standardization
    % Eq. 21: Jl*^T{diag(S_J)}^-1*Jl
    JJ1 = bsxfun ( @times, J1.^2, 1./diag(S_j));
%     for i = 1:Nd
%       JJ1(:,i) = (1/S_j(i,i))*J1(i,:).^2;
%     end
  elseif isequal(alpha,'timevarying')
    [U,s,~] = csvd(G); % standard-form regularization
    for k = 1:Nt
      alpha1 = l_curve(U,s,data(:,k),'tikh');
      
      % standardiaztion
      T = K'*pinv(K*K' + alpha1*H); % Eq. 6 or Eq.11: T = K^TH[HKK^T + alpha*H]+
      S_j = T * K;  % Eq. 17: S_j = K^TH[HKK^T + alpha*H]+
      J1 = S_j * (T * phi(:,k));
      
      % Using Dale et al.standardization
      % Eq. 21: Jl*^T{diag(S_J)}^-1*Jl
      JJ1 = bsxfun ( @times, J1.^2, 1./diag(S_j));
%       for i = 1:Nd
%         JJ1(k,i) = (1/S_j(i,i))*J1(i,:).^2;
%       end
    end
  else
    error('alpha is invalid...');
  end
elseif isnumeric(alpha)
  % standardiaztion
  T = K'*pinv(K*K' + alpha*H); % Eq. 6 or Eq.11: T = K^TH[HKK^T + alpha*H]+
  S_j = T * K;  % Eq. 17: S_j = K^TH[HKK^T + alpha*H]+;
  J1 = S_j * (T * phi);
  
  % Using Dale et al.standardization
  % Eq. 21: Jl*^T{diag(S_J)}^-1*Jl
  JJ1 = bsxfun ( @times, J1.^2, 1./diag(S_j));
%   for i = 1:Nd
%     JJ1(:,i) = (1/S_j(i,i))*J1(i,:).^2;
%   end

end

if flag1 == 1
  JJ2 = JJ1(1:Nd/3,:)+JJ1((Nd/3)+1:(2*Nd/3),:)+JJ1((2*Nd/3)+1:Nd,:);
else
  JJ2 = JJ1;
end

% -----------------------------------
% Eq. 20: Jl*^T{S_J}^-1*Jl
% JJ2 = zeros(Nt,Nd/3);
% s = 0;
% for i=1:3:Nd  %
%     s = s+1;
%     inv1(:,:,s)= inv(diag(diag(S_j(i:3*s,i:3*s))));
%     JJ2(:,s)= diag(J1(i:3*s,:)'*inv1(:,:,s)*J1(i:3*s,:));
% end
% ------------------------------------
% test only
% JJ2 = zeros(Nt,Nd/3);
% s = 0;
% for j = 1:3:Nd
%   s = s+1;
%   JJ2(:,s) = JJ1(:,j) + JJ1(:,j+1) + JJ1(:,j+2);
% end
