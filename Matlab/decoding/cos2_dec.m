%function [s,r,m]=cos2_dec(ch,ssc,fpos,cwidth,mfl)
%
%Decode signal values in cos^2 channel representation.
%
%Each column of CH is taken to be a channel vector,
%and all columns of CH are decoded in parallel.
%Corresponding columns in S and R, contain lists of
%signal values and relevance measures, both sorted
%according to the relevance, R.
%
%  CH     Channel value matrix [CHANNELS x SAMPLES]
%         or column vector if just one sample.
%  SSC    channel distance (default 1)
%  FPOS   first channel position (default 1)
%  CWIDTH Channel width in terms of cos^2 frequency
%         (default pi/3)
%  MFL    modular flag (default 0)
%
% Returns
% S    Decoded positions
% R    Channel sum relevance
% M    Vector magnitude relevance
%
%NOTE: Multiple signal values (s) are always
%      returned. Signal values with zero
%      relevance (r) should be ignored.
%
%See also COS2_ENC and COS2_POS
%
%Per-Erik Forssen, 2001

function [S,R,M]=cos2_dec(ch,ssc,fpos,cwidth,mfl)

if nargin<2 | isempty(ssc),ssc=1;end;
if nargin<3 | isempty(fpos),fpos=1;end;
if nargin<4 | isempty(cwidth),cwidth=pi/3;end;
if nargin<5 | isempty(mfl),mfl=0;end;

ocnt=pi/cwidth-1;
if ((ocnt-floor(ocnt))~=0)|(ocnt<1),
  error('cwidth has to be pi/N where N is an integer >= 2');
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
S=S(ind+repmat(([1:samples]-1)*rows,rows,1));
M=M(ind+repmat(([1:samples]-1)*rows,rows,1));
