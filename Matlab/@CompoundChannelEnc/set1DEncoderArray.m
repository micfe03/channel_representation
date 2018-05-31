%CompoundChannelEnc/set1DEncoderArray
%
%  obj = set1DEncoderArray(obj, enc1Ds)
% 
%  Gives the user complete control over the CompoundChannelEnc class by
%  submitting a custom enc1Ds cell array. This is an array of ChannelEnc1D
%  (or compatible) objects and has the same structure as 'chanlayout', i.e. 
%  each row represents one layer, each column represents an input
%  dimension, an empty cell means "don't care", channels from the same row
%  are outer-product-ed, different rows are concatenated.
%
%  [CVL 2007 (see README.txt for credits)]

