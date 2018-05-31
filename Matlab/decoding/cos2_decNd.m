function [s,r] = cos2_decNd(ch,ssc,fpos,cwidth,mfl,Nvalues)

% function [s,r] = cos2_decNd(ch,ssc,fpos,cwidth,mfl,Nvalues)
%
% Channel decoding of an N-dimensional channel representation.
%
%   s       N/Nvalues matrix of decoded values
%   r       1/Nvalues matrix of certainties
%
%   ch      M1/M2/.../MN Channel value array, where Mn is the number
%           of channels in the n:th dimension
%   ssc     N/1 channel distances (default 1)
%   fpos    N/1 first channel positions (default 1)
%   cwidth  N/1 Channel widths in terms of cos^2 frequency
%           (default pi/3)
%   mfl     N/1 modular flags (default 0)
%   Nvalues Number of decoded values (default 1)
%
%
%
% bjojo, 2004-05-28

sz = size(ch);
Ndim = length(sz);

if nargin<2 | isempty(ssc)
  ssc = ones(Ndim,1);
end
if nargin<3 | isempty(fpos)
  fpos = ones(Ndim,1);
end
if nargin<4 | isempty(cwidth)
  cwidth = pi/3*ones(Ndim,1);
end
if nargin<5 | isempty(mfl)
  mfl = zeros(Ndim,1);
end
if nargin<6 | isempty(Nvalues)
  Nvalues = 1;
end

% Special case
if Ndim==2 & (sz(1)==1 | sz(2)==1)
  if sz(1)==1
    ch = ch';
  end
  [s r] = cos2_dec(ch,ssc,fpos,cwidth,mfl);
  s = s(1:Nvalues);
  r = r(1:Nvalues);
  return
end

% Sum over local regions
%========================

chw = comp_cert(ch,Ndim,cwidth,mfl);

% Decode local maxima regions
%=============================

nms = non_max_suppressionNd(chw,mfl);
[ind foo r] = find(nms(:));
[r I] = sort(r,1,'descend');
ind = ind(I);

if Nvalues>length(r)
  Nvalues = length(r);
end
r = r(1:Nvalues);

s = zeros(Ndim,Nvalues);
psz = [1 cumprod(sz)];
for nv=1:Nvalues
  sub = my_ind2sub(sz,ind(nv));

  % Channel decode
  %================
  
  for n=1:Ndim
    N = round(pi/cwidth(n));
    ind_n = (1:N)-ceil(N/2);
    if mfl(n)==0
      indind = ((sub(n)+ind_n)>0) & ((sub(n)+ind_n)<=sz(n));
      ind_n = ind_n(indind);
      IND_n = ind(nv) + psz(n)*ind_n;
      ch_n = zeros(N,1);
      ch_n(indind) = ch(IND_n);      
    else    
      ind_n = mod(ind_n+sub(n)-1,sz(n))+1-sub(n);
      IND_n = ind(nv) + psz(n)*ind_n;
      ch_n = ch(IND_n)';
    end
    
    fpos_n = fpos(n) + (sub(n)-2)*ssc(n);
    
    sn = cos2_dec(ch_n,ssc(n),fpos_n,cwidth(n),0);
    s(n,nv) = sn(1);
    if mfl(n)
      s(n,nv) = mod(s(n,nv)-fpos(n),ssc(n)*sz(n))+fpos(n);
    end
  end
end

%if Ndim==2 & sz(2)==1
%  s = s(1);
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IND = my_ind2sub(siz,ndx)

% Same as ind2sub.m except that the output is a vector.
%
% IND - N/M matrix of indeces
%
% siz - 1/N vector of size in each dimension
% ndx - 1/M vector of indeces

n = length(siz);
k = [1 cumprod(siz(1:end-1))];
ndx = ndx - 1;
IND = zeros(length(siz),length(ndx));
for i = n:-1:1,
  IND(i,:) = floor(ndx/k(i))+1;
  ndx = rem(ndx,k(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function chw = comp_cert(ch,Ndim,cwidth,mfl)

chw = ch;
if nnz(round(pi./cwidth)==3)==Ndim
  % Faster code for case cwidth=pi/3
  % Sum over 3x3x...x3-regions
  w = ones([3 ones(1,Ndim-1)]);
  for n=1:Ndim
    if mfl(n)==0
      chw = imfilter(chw,w,'same','corr');
    else
      chw = imfilter(chw,w,'same','corr','circular');
    end  
    w = shiftdim(w,-1);
  end
else
  % Slower code for general case of cwidth
  % Sum over N1xN2x...xNm-regions
  for n=1:Ndim
    M = round(pi/cwidth(n));
    w = ones([M ones(1,Ndim-1)]);
    for k=1:n-1
      w = shiftdim(w,-1);
    end
    
    if mfl(n)==0
      chw = imfilter(chw,w,'same','corr');
    else
      chw = imfilter(chw,w,'same','corr','circular');
    end
  end 
end
