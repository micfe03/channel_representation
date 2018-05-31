%CompoundChannelEnc/encodeImg
%
%  chimg = encodeImg(enc, img)
%
%  Similar to ChannelEnc1D/encode, but expects 3D array 'img', where each 
%  pixel position contains 'N' values to be encoded. Reshapes the image,
%  encodes each pixel, and reshapes back.
%
%  [CVL 2007]
