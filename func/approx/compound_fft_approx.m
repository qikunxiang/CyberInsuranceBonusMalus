function p_cp = compound_fft_approx(p, freq_pgf, tilt)
%COMPOUND_FFT_APPROX Approximate the pmf of a discrete compound
%distribution using fast Fourier transform
% Inputs: 
%       p: the pmf of the severity distribution
%       freq_pgf: the pgf of the frequency distribution
%       tilt: the tilting parameter used for improving the accuracy of
%           approximation, default is 0
% Outputs:
%       p_cp: the pmf of the compound distribution

if ~exist('tilt', 'var') || isempty(tilt)
    tilt = 0;
end

atoms = (1:size(p, 1))';
p_cp = ifft(freq_pgf(fft(p .* exp(-tilt * atoms)))) .* exp(tilt * atoms);
p_cp = real(p_cp);

end

