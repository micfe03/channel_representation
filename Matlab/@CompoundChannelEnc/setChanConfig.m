%CompoundChannelEnc/setChanConfig
%
%  obj = setChanConfig(enc1D, mode, chanlayout, bounds, cscales, mflags)
%
%  Sets the number of channels, width and position in input space. Some
%  different modes are available for defining the channel positions
%
%  Parameters: 
%    enc1D:      Template for the 1D encoders. Must be equipped with the
%                desired basis function.
%
%    mode:       'exterior' or 'interior' channel placement strategy (see
%                ChannelEnc1D/setChanConfig)
%
%    chanlayout: The "channel layout matrix". Each column corresponds 
%                to one dimension of the data, and each row specifies
%                one "encoding layer". See ChannelEncND for furher info.
%
%   bounds:      Nx2 matrix with bounds for each input dimension.
%
%   cscales:     N-vector with channel scaling for each input dimension
%                (optional, default: all ones)
%
%   mflags:      N-vector with modular flag for each input dimension.
%                (optional, default: all zeros)
%
%   If some parameter is [], it uses the default value.
%
%   Example:
%   > nchans = [8 8 8]; 
%   > bounds = [0 1 ; 0 2*pi ; 0 1];  
%   > mflags = [0 1 0];
%   > enc = setChanConfig(enc1D, 'exterior', nchans, bounds, [], mflags);
%
%  [Erik Jonsson, 2007]
