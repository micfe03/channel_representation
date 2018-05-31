function [S, R, M] = cos2decode(ch, cscale, mfl, nmodes)
%
%  [S, R, M] = cos2decode(channel, cscale, mflag, nmodes)
%
%  Standard decoding of cos^2 kernel (cos2kernel.m). This function follows 
%  the standard channel decoding interface used by ChannelEnc1D 
%  (see ChannelEnc1D/decde)
%
%  Parameters:
%    ch:       Matrix of channel vectors as returned from ChannelEnc1D/encode.
%              Each column is one channel-encoded value
%
%    cscale:   Channel scaling. Default is 1, corresponding to kernel 
%              cos(pi/3*x)^2.
%
%    mflag:    True if the domain is modular
%
%    nmodes:   Number of modes to return. Optional, default is 1. 
%              (not implemented yet)
%
%  Return values:
%    S:        Decoded signal values
%    R:        Certainty (channel sum relevance)
%    M:        Vector magnitude relevance (not required/returned by the 
%              ChannelEnc1D/decode interface, but present in Per-Erik's
%              original code.
%
%  [CVL 2007 (see README.txt for proper credits)]


% Using Per-Erik's old code, but sets ssc and fpos to constants. This way,
% the code looks as similar as possible to the old cos2_enc.
ssc = 1;
fpos = 1;

% Argument handling
if nargin<2 || isempty(cscale), cwidth=pi/3; end;
if nargin<3 || isempty(mfl), mfl=0; end;
if nargin<4 || isempty(nmodes), nmodes=1; end;

% Per-Erik's decoding code uses cwidth, which is in radians, such that
% the kernel is cos(cwidth*x)^2. However, now we want to be more general
% and instead specify a scaling related to some standard kernel. To keep 
% the rest of the code as similar as possible, we here create the old variable 
% cwidth from cscale.
cwidth = pi/3/cscale;

ocnt=pi/cwidth-1;
if ((ocnt-floor(ocnt))~=0)|(ocnt<1),
  error('cscale must be N/3, where N is an integer >= 2');
end;

chans=size(ch,1);
samples=size(ch,2);
%
% We now pad the channel vector to be able to loop directly
% over it. For modular channel sets we copy the channel values
% from the other side of the vector, otherwise we pad with zeros.
%
if mfl,
  ch2=[ch(chans-ocnt+2:chans,:);ch;ch(1:ocnt-1,:)];
else
  ch2=[zeros(ocnt-1,samples);ch;zeros(ocnt-1,samples)];
end;

%
% Local inverse for [ch(k) ch(k+1) ... ch(k+w) ]
%    Solutions outside [w/2-0.5,w/2+0.5] are better
%    described by another group of coefficients
%
llim=ocnt/2-.5;
ulim=ocnt/2+.5;
S=zeros(chans,samples);
R=zeros(chans,samples);
M=zeros(chans,samples);
if(ocnt==1),
  % The not very useful case of cwidth=pi/2
  for k=1:chans-1,
    offset=atan(sqrt(ch2(k+1,:)./(ch2(k,:)+eps)))*2/pi;
    S(k,:)=k+offset;
    if(chans>2)      
      switch(k)
       case 1,
	valfl=(ch2(k,:)>=ch2(k+2,:));
       case chans-1,
	valfl=(ch2(k+1,:)>=ch2(k-1,:));
       otherwise
	valfl=(ch2(k,:)>=ch2(k+2,:)).*(ch2(k+1,:)>=ch2(k-1,:));
      end
    else
      valfl=ones(2,samples);
    end
    R(k,:)=(ch2(k,:)+ch2(k+1,:)).*valfl;
    M(k,:)=abs(ch2(k,:)+i*ch2(k+1,:)).*valfl;
  end
else
  a=2*cwidth*[0:ocnt]';
  W=2*pinv([cos(a) sin(a) ones(ocnt+1,1)]);

  for k=1:chans+ocnt-2,
    cvec=W*ch2(k:k+ocnt,:);
    offset=atan2(cvec(2,:),cvec(1,:))/2/cwidth;
    S(k,:)=k-ocnt+1+offset;
    valfl=(offset>llim).*(offset<=ulim);
    R(k,:)=cvec(3,:).*valfl;
    M(k,:)=abs(cvec(1,:)+i*cvec(2,:)).*valfl;
  end;
end

%
% Add scaling and offset
%
% FIXME: Borde man helt enkelt kunna ta bort?
if mfl,
  S=mod(S-1,chans)*ssc+fpos;
else
  S=(S-1)*ssc+fpos;
end
S=S.*(R>0);                         % Zero value for zero relevance

%
% Sort according to relevance
%

%Dominant scalar only
%[R,ind]=max(R,[],1);
%S=S(ind+([1:samples]-1)*chans);

[R,ind]=sort(-R);
R=-R;
rows=size(ind,1);

S = S(ind+repmat(([1:samples]-1)*rows,rows,1));
M = M(ind+repmat(([1:samples]-1)*rows,rows,1));

% Return only the 'nmodes' strongest modes
R = R(1:nmodes, :);
S = S(1:nmodes, :);
M = M(1:nmodes, :);


