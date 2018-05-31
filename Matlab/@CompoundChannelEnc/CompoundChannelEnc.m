classdef CompoundChannelEnc
%
%  CompoundChannelEnc
%
%  Implements multi-dimensional channel encoding by combining
%  one-dimensional encodings using any combination of outer products and
%  concatenation. 
%
%  The key to this is the 'chanlayout' matrix. Each row in this matrix
%  specifies one 'layer', and defines how many channels to use in each
%  dimension. Each dimension is encoded separately, and the outer product 
%  is taken between all these encodings. A zero means 'ignore this
%  dimension'. Finally, all layers are concatenated. This is perhaps best
%  illustrated by some examples:
%
%  Ex1: chanlayout = [8 8 8];
%    Encodes a 3D vector with the outer product of 8 channels in each 
%    dimension. Produces 8^3=512 channel coefficients.
%         
%  Ex2: chanlayout = 8*eye(3);
%    Encodes a 3D variable with 8 channels in each dimension separately and 
%    concatenates the results. Produces 3*8=24 coefficients
%        
%  Ex3: chanlayout = [5 5 5 ; 10*eye(3)] 
%    First takes the outer product of 5 channels in each dimension. Then
%    concatenates this with 10 channels separate for each dimension.
%    Produces 5^3+3*10 = 155 coefficients.
%
%  The call to the constructor defines the basis function, similarily to
%  ChannelEnc1D. However, it is NOT possible (yet?) to include the channel
%  placement parameters to the constructor. Instead, call setChanConfig
%  directly after the constructor call (type 'help CompoundChannelEnc/setChanConfig'
%  for more info.
%
%  Internally, the multi-dimensional encoder is represented as a number of
%  encoders stored in a cell array similar to the 'chanlayout' array. This
%  array can be set using set1DEncoderArray
%
%  Methods:
%    setChanConfig
%    encode
%    encodeImg
%    set1DEncoderArray
%
%
%  Ex: Create a 8*8*8 outer-product encoding with B2-kernels and one
%      modular domain
%  > enc = CompoundChannelEnc('bsp2'); 
%  > enc = setChanConfig(enc, 'exterior', [8 8 8], [0 1 ; 0 2*pi ; 0 1], [0 1 0]);
%  > ch = encode(enc, [0.2;pi;0.1]);
%
%  [CVL 2007 (see README.txt for credits)]


%  Each of these is responsible for encoding a scalar value, using it's
%  own settings for basis function, bounds, number of channels etc.
%
%  It would be possible to add some interface to this class that lets
%  a user specify this cell array completely. 

  
properties (SetAccess = 'private', GetAccess = 'public')
  
  enc_template; % Template 1D encoder (created with the correct kernel)

  enc1Ds;       % Cell array of one-dimensional encoders, with structure
                % corresponding to chanlayout

  decf;         % Decoder (not supported yet)
  
  % Channel encoding parameters 
  chanlayout = []; % Number of channels in each layer
  nchans = [];     % Total number of channels
  bounds = [];
  mflags = [];
  cscales = [];
  
end
  
  
methods
  
  % Constructor
  function obj = CompoundChannelEnc(varargin)
    obj.enc_template = ChannelEnc1D(varargin{:});
  end

  
  
  function obj = setChanConfig(obj, mode, chanlayout, bounds, cscales, mflags)
    
    % Default input arguments
    if nargin<5, cscales = ones(1,size(bounds,1)); end;
    if nargin<6, mflags = zeros(1,size(bounds,1)); end;

    obj.enc1Ds = cell(size(chanlayout));
    obj.chanlayout = chanlayout;
    obj.bounds = bounds;
    obj.cscales = cscales;
    obj.mflags = mflags;
    
    for ri = 1:size(obj.enc1Ds,1)
      for ci = 1:size(obj.enc1Ds,2)

        if chanlayout(ri,ci)~=0
          % Create a new copy of the template encoder with different
          % parameters
          obj.enc1Ds{ri,ci} = setChanConfig(obj.enc_template, mode, ...
                chanlayout(ri,ci), bounds(ci,:), cscales(ci), mflags(ci));
        end
        
      end
    end

    dummy = obj.chanlayout;
    dummy(dummy==0) = 1;
    obj.nchans = sum(prod(dummy,2));
    
  end

  
  function obj = set1DEncoderArray(obj, enc1Ds)
    obj.enc1Ds = enc1Ds;

    % Make sure that the other fields are consistent
    obj.chanlayout = [];
    obj.bounds = [];  % Not applicable - may be different in each layer
    obj.enc_template = []; % Not applicable
    obj.cscales = []; % Not applicable - may be different in each layer
    obj.mflags = [];  % Not applicable - may be different in each layer 
                      % (even if that would be stupid)

    for ri = 1:size(obj.enc1Ds,1)
      for ci = 1:size(obj.enc1Ds,2)

        if ~isempty(enc1Ds{ri,ci})
          obj.chanlayout(ri,ci) = enc1Ds{ri,ci}.nchans;
        else
          obj.chanlayout = 0;
        end
        
      end
    end

    dummy = obj.chanlayout;
    dummy(dummy==0) = 1;
    obj.nchans = sum(prod(dummy,2));
  end

  

  % --- Encoding ---------------------------------------------------------
  function chtot = encode(obj, vals)
  
    % Some error checking
    if size(vals,1)~=size(obj.enc1Ds,2);
      error('Dimensionality mismatch in encode');
    end
    
    % Main loop
    chtot = [];
    for ri = 1:size(obj.enc1Ds, 1) % Loop through the layers
      
      ch1 = ones(1, size(vals, 2));
      for di = 1:size(obj.enc1Ds, 2)  % Loop through the signal dimensions
        
        if not(isempty(obj.enc1Ds{ri, di}))
          % Encode using encoder (ri,di) and form column-wise outer 
          % product with previous
          ch = encode(obj.enc1Ds{ri,di}, vals(di,:));
          ch1 = kroncomp(ch, ch1);
        end
      end
      
      % Gather in a matrix
      chtot = [chtot ; ch1];
    end
    
  end

  
  function chimg = encodeImg(obj, img)
    % Rehape the image to a matrix, where each column is from one pixel.
    % Encode each pixel, and reshape back

    sz = size(img);
    
    img = reshape(img, sz(1)*sz(2), []);
    chimg = encode(obj, img');
    chimg = reshape(chimg', sz(1), sz(2), []);
  
  end

  
  function vals = decode(obj, ch)
    error(['Decoding not implemented. Type help CompoundChannelEnc/decode for more info']);
  end
  
  function vals = decodeImg(obj, ch)
    error(['Decoding not implemented. Type help CompoundChannelEnc/decode for more info']);
  end

   
end

end
