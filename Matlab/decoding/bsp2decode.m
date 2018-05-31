function [S, E, R] = bsp2decode(channel, cscale, mflag, nmodes)
%
%  [vals, cert] = bsp2decode(channel, cscale, mflag, nmodes)
%
%  Standard decoding of second order B-splines channel encoding
%  (bsp2kernel). This function follows the standard channel decoding
%  interface used by ChannelEnc1D (see ChannelEnc1D/decde)
%
%  Parameters:
%    channel:  Matrix of channel vectors as returned from ChannelEnc1D/encode.
%              Each column is one channel-encoded value
%
%    mflag:    True if the domain is modular
%
%    nmodes:   Number of modes to return. Optional, default is 1. 
%
%  [CVL 2007-2012 (see README.txt for proper credits)]

% Default input arguments
if nargin<2, mflag=0; end;
if nargin<3, nmodes=1; end;

% This function is somewhat restricted from cos2enc. Check for unsupported
% arguments (please fix in the future!)
if cscale~=1
  error('bsp2decode only supports scale=1');
end
  
K=size(channel,1);   % Number of channels
s1d=size(channel,2); % Number of values

% channel = channel./(ones(s1d,1)*sum(channel,1)); % this is a fix that is required to get the absolute level correct

if mflag
  % Periodic decoding, from bscd_per3 by MFE
  
  R=sum(channel);
  
  % iir filter for spline interpolation
  rv=zeros(K,1);
  rv(1)=3/4;
  rv(2)=1/8;
  rv(K)=1/8;
  RV=fft(rv);
  C=fft(channel,[],1);
  C=(RV.^-1*ones(1,s1d)).*C;
  channel=ifft(C,[],1);
  channel=[channel; channel(1:4,:)];
  
  % zeros of influence function
  n=conv2(channel,[1 -2 0 2 -1]','valid');
  p=conv2(channel,[-1 0 2 0 -1]','valid')./2;
  q=conv2(channel,[1 6 0 -6 -1]','valid')./4;
  S=p.^2-q.*n;
  S=(p-sqrt((S>0).*S))./(n+eps);
  msk=abs(S)<=0.5;
  
  % energy at zeros
  E=conv2(channel,[1/48 1/2 23/24 1/2 1/48]','valid')+S.*(q+S.*(-p+S.*n./3))./2;
  %E=23/24-channel(3:K,:).*(B3spline(b1-.5)+B3spline(b1+.5))-channel(2:K-1,:).*(B3spline(b1+.5)+B3spline(b1+1.5))...
  %-channel(1:K-2,:).*B3spline(b1+1.5)-channel(4:K+1,:).*(B3spline(b1-1.5)+B3spline(b1-.5))-channel(5:K+2,:).*B3spline(b1-1.5);
  E(~msk)=0;
  
  if nmodes==1
      % strongest mode selection
      [E,ind]=min(23/24-E);
      %S1=S(ind+[0:s1d-1]*K)+ind;
      S=S(ind+[0:s1d-1]*K);
      S=mod(S+ind+2,K);
  else
      [E,ind]=sort(23/24-E);
      E=E(1:nmodes,:);
      S=S(ind+repmat([0:s1d-1]*K,K,1));
      S=mod(S+ind+2,K);
      S=S(1:nmodes,:);
  end
  %R=R-24/23*E;

else
  % Linear decoding, adapted from bscd_lin v3.3 by MFE, PEF
    
  % iir filter for spline interpolation
  z=2*sqrt(2)-3;
  channel=filter([1],[1,-z],[zeros(1,s1d); channel; zeros(1,s1d)],0);
  zi=z^3/(z^2-1)*channel(K+2);
  channel=filter([-z],[1,-z],channel(K+2:-1:1,:),zi);
  channel=8*channel(K+2:-1:1,:);
  
  % zeros of influence function
  n=conv2(channel,[1 -2 0 2 -1]','valid');
  p=conv2(channel,[-1 0 2 0 -1]','valid')./2;
  q=conv2(channel,[1 6 0 -6 -1]','valid')./4;
  S=p.^2-q.*n;
  S=(p-sqrt((S>0).*S))./(n+eps);
  msk=abs(S)<0.505; % hack: 1% is allowed
  
  % energy at zeros
  E=23/24-conv2(channel,[1/48 1/2 23/24 1/2 1/48]','valid')-S.*(q+S.*(-p+S.*n./3))./2;
  %E=23/24-channel(3:K,:).*(B3spline(b1-.5)+B3spline(b1+.5))-channel(2:K-1,:).*(B3spline(b1+.5)+B3spline(b1+1.5))...
  %-channel(1:K-2,:).*B3spline(b1+1.5)-channel(4:K+1,:).*(B3spline(b1-1.5)+B3spline(b1-.5))-channel(5:K+2,:).*B3spline(b1-1.5);
  E(~msk)=23/24;
  
  if size(E,1)>1
      if nmodes==1
          [E,ind] = min(E);
      else
          [E,ind]=sort(E);
          E=E(1:nmodes,:);
      end
  else
      ind = 1;
  end
      
  S=S.*msk; % wrong offsets are replaced with channel center
  if nmodes==1
      S=S(ind+(0:s1d-1)*(K-2))+ind+1;
  else
      S=S(ind+repmat([0:s1d-1]*(K-2),K-2,1))+ind+1;
      S=S(1:nmodes,:);
  end
  
  %R = 'dummy'; % FIXME: Just to be able to compare and see what happens
  
end

E = 23/24 - E; % fix by mfe: return confidence, not error!


% Reshaping from image to matrix first
%K=size(channel,3);
%2d=size(channel(:,:,1));
%s1d=prod(s2d);
%channel=reshape(channel(:),[s1d K])';


%  if 0
%    % strongest mode selection
%    if nargin<3
%      [E,ind]=min(E);
%    else
%      E=E(ind+(0:s1d-1)*(K-2));
%    end
%  end

