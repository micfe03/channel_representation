function [S, R, M] = rectdecode(ch, cscale, mfl, nmodes)
%
%  [S, R, M] = rectdecode(channel, cscale, mflag, nmodes)
%
%  Standard decoding of rect kernel (tectkernel.m). This function follows 
%  the standard channel decoding interface used by ChannelEnc1D 
%  (see ChannelEnc1D/decode)
%
%  Parameters:
%    ch:       Matrix of channel vectors as returned from ChannelEnc1D/encode.
%              Each column is one channel-encoded value
%
%    cscale:   Channel scaling. Default is 1, corresponding to kernel 
%              rect(x).
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

[R S] = max(ch);


