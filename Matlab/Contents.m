%% CVL Channel Representation Toolbox
%
% File revision $Id: Contents.m 3638 2006-06-20 13:48:13Z jowi $
%
% This toolbox contains functions for channel encoding and decoding, using
% various basis functions, e.g. bspline, cos2, p-channels etc. Learning
% methods based on channel representations do not belong here, but in
% pattern_recog.
%
% Previously, there were many different channel encoders / decoders with
% slightly different syntax and semantics. In a unification attempt, these
% have been (partly) replaced by two channel encoder classes. The old
% functions are still available under 'obsolete'
%
% The classes can use any of the kernels and decodings available and should
% be the primary interface to the toolbox. These classes then use different
% kernels and decoders internally. However, the kernels and decoders are
% also accessible from the outside.
%
% The idea behind the object orientated structure is that you can create
% your channel encoders once and for all, in the top of your scripts, and
% then just use 'encode' and 'decode' etc, regardless of the details about
% your encoding.
%
% Note that the classes are written using the new, unofficial Matlab OO 
% system (MCOS), and requires Matlab 7.3. The classes will be maintained to 
% work with the latest matlab version until MCOS is officially released.
%
% *Classes:*
%
%   ChannelEnc1D            - Class for encoding/decoding of scalars
%
%   CompoundChannelEnc      - Class for computing compound encodings using
%                             outer-products and concatenation of 1D encodings
%
% *Kernel functions:*
%
%    bkspline               - N'th order B-spline kernel
%    bsp2kernel             - Second order B-spline kernel
%    bsp3kernel             - Third order B-spline kernel
%    cos2kernel             - cos^2 kernel
%    pkernel                - P-channel kernel(s)
%    rectkernel             - Regular histogram kernel (box function)
%    bsp2kernel_diff        - The derivative of the B2-kernel
%      (more derivatives of kernels can be added)
%
% *Decoders:*
%
%    bsp2decode             - Decoding of B2-encodings
%    cos2decode             - Decoding of cos2-encodings
%    pdecode                - Decoding of P-channel encodings
%
% *Optimized spatio-featural representations*
% 
%    pwpolyenc_mex          - Spatio-featural representation by piecewise
%                             polynomials (P-channels as special case)
%
%    pwpoly2bspline         - Converts piecewise polynomial representations
%                             to b-spline representations
%

%%
cvlUseToolbox('testimages');  % For baboon

%% ChannelEnc1D
help ChannelEnc1D/ChannelEnc1D

%% ChannelEnc1D/setChanConfig
help ChannelEnc1D/setChanConfig

%% ChannelEnc1D/encode
help ChannelEnc1D/encode

%% ChannelEnc1D/decode
help ChannelEnc1D/decode

%% ChannelEnc1D/centers
help ChannelEnc1D/centers

%% ChannelEnc1D/basisMatrix
help ChannelEnc1D/basisMatrix

%% ChannelEnc1D/encodeImg
help ChannelEnc1D/encodeImg

%% ChannelEnc1D/decodeImg
help ChannelEnc1D/decodeImg

%% ChannelEnc1D/encodeDensity
help ChannelEnc1D/encodeDensity

%% --ChannelEnc1D Examples--
enc = ChannelEnc1D('cos2', 'exterior', 6, [0 1]);
ch = encode(enc, [0 0.3 1])  % Encode 3 values
vals = decode(enc, ch)
c = centers(enc)


%%
figure(1);
plot(basisMatrix(enc, 1000))  % Plot the basis functions

figure(2);
img = im2double(rgb2gray(imread('baboon.tif')));
chimg = encodeImg(enc, img);           % Encode an entire image
imagesc(chimg(:,:,3)); colormap gray;  % Display channel 3

figure(3);
img2 = decodeImg(enc, chimg);
imagesc(img2); colormap gray;  % Display the decoded image



%% CompoundChannelEnc
help CompoundChannelEnc/CompoundChannelEnc

%% CompoundChannelEnc/setChanConfig
help CompoundChannelEnc/setChanConfig

%% CompoundChannelEnc/set1DEncoderArray
help CompoundChannelEnc/set1DEncoderArray

%% CompoundChannelEnc/encode
help CompoundChannelEnc/encode

%% CompoundChannelEnc/decode
help CompoundChannelEnc/decode

%% CompoundChannelEnc/encodeImg
help CompoundChannelEnc/encodeImg

%% --CompoundChannelEnc examples--

% Concatenation of 5 channels in each dimension
enc = CompoundChannelEnc('bsp2');
enc = setChanConfig(enc, 'exterior', 5*eye(3), [0 1 ; 0 1 ; 0 1]); 
ch = encode(enc, [0 0.3 ; 0 0 ; 0 1])  % Encode 2 3D-values


%% bkspline
help bkspline

%% bsp2kernel
help bsp2kernel             

%% bsp3kernel
help bsp3kernel

%% cos2kernel
help cos2kernel

%% pkernel
help pkernel

%% rectkernel
help rectkernel

%% bsp2kernel_diff
help bsp2kernel_diff


%% pwpolyenc_mex
% help pwpolyenc_mex (coming soon)

%% pwpoly2bspline
% help pwpoly2bspline (coming soon)


