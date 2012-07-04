% Rice Wavelet Toolbox
% Version 2.4   Dec 2002
%
% Wavelet Transforms.
%     mdwt - Discrete orthogonal wavelet transform using Mallat alg. (1D and 2D)
%     midwt - Inverse discrete orthogonal wavelet transform
%     mrdwt - Undecimated (redundant) discrete wavelet transform (1D and 2D)
%     mirdwt - Inverse undecimated discrete wavelet transform
%     daubcqf - Daubecheis filter coefficients
%
% Wavelet Domain Processing.
%     denoise - Denoise signal and images by thresholding wavelet coefficients
%     HardTh - Hard thresholding
%     SoftTh - Soft thresholding
%
% Other.
%     makesig - Create Donoho-Johnstone test signals
%     compile - compile the Rice Wavelet toolbox

% Added by Robert Brockman II, April 2012:
%     test_*dwt.m - Unit Tests for above wavelet transforms.  Requires
%                   xUnit in path, use "runtests" in MATLAB Command Window to access.
%     lena512.mat - Standard "Lena" grayscale test image.
