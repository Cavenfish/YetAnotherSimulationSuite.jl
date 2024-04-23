
The best way to explain my VACF code/procedure is to explain how each parameter influences the results. So first I will list what parameters are available to vary.

1) Normalization of auto-correlation
2) Window to be used
3) Degree of zero-padding
4) Mirroring the ACF

## Normalization

## Window

The window function is extremely important, and unfortunately there is no one size fits all window. Some are know to generally work well for all cases, like the Hann window. However, its best to always double check that the window in use does not cause spectral leakage.  Below I give some windows available, with the convention that $n$ is the index and $N$ is the total length.

#### Hann
$$
W_\text{Hann}(n) = \sin^2(\frac{n\pi}{N})
$$
#### Hann Mirroring Version
$$
W_\text{HannM}(n) = \cos^2(\frac{n\pi}{2(N-1)})
$$
#### Welch
$$
W_\text{Welch} = 1 - (\frac{n - \frac{N}{2}}{\frac{N}{2}})^2
$$

## Zero-padding

This is just adding zeros to end of the ACF array. It is an artificial way to smooth the FFT of the ACF. Keep in mind that this smoothing does not actually give new information, rather it simply fits a sin wave between points of data you had previously. Any gain in resolution from the zero-padding is essentially just interpolation points and not real resolution. 

## Mirroring the ACF

This is similar to zero-padding, it will attempt to smooth the FFT of the ACF but will only be an interpolation. This comes with an added increase in smoothness (and caveat), since you mirror your ACF after zero-padding it will double the interpolation resolution of your zero-padding. The caveat here is that you need to ensure your window function is made for mirroring (ie. it starts at 1 rather than $\approx$ 0). This is important since the signal you want to perform an FFT on needs to go from 0 to 1 to 0, and if you mirror your ACF then you need to start from 1 to get that signal format. Lastly, its easy to not make this mistake but its costly. The mirrored ACF cannot have a repeated value within it, this mistake usually happens at the joining of the mirror (ie. a repeated 1 value). This will cause spectral artifacts that are more pronounced with larger time steps.  