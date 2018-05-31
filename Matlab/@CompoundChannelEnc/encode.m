%CompoundChannelEnc/encode
%
%  ch = encode(obj, vals)
%
%  Encodes all values in 'vals' in multiple subspaces using the combination
%  of 1D channel encodings specified in 'obj'. 
%
%  Parameters:
%    obj:   ChannelEncND object
%    vals:  Matrix of vectors to encode. Each column is one vector
%
%  Each column in the returned matrix 'ch' is the channel encoding of each
%  column of the input. The first dimension in the input corresponds to
%  the least significant part of the serialized index. We can think of 
%  it like:
%    > ch = encode(enc, [aval bval cval])
%    > ch = reshape(ch, [na, nb, nc]);
%    > ch(ai, bi, ci)
%  i.e, think about 'ch' as a serialized 3D-array of dimensions (na,nb,nc)
%
%  [Erik Jonsson, 2007]
