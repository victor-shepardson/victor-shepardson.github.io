Note: this webpage needs to fetch audio from soundcloud.com, and scripts from cdn.mathjax.org for math rendering

# Data Driven Gestural Control of a Synthesis Process by an Audio Stream

## Project Milestone

### Victor Shepardson | Dartmouth College | CS 174 Fall 14

## Problem

Given a synthesis process $S$ which converts a time series of synth parameters to a time series of audio samples, a basis function $b$ which maps audio to a perceptual space, and an input audio series $a$, we attempt to find the parameter series $y$ which minimizes $\lvert b(a) - b(S(y)) \rvert$. The result is a kind of gestural control of the synthesis process by the audio stream.

## Motivation

The intended use of this system is artistic. The idea is to let a performer or composer map the essential gestures of a performance into a different timbral space. Exactly what constitutes a gesture (pitch? spectral envelope?) depends on $b$ and what timbres can be produced depend on $S$. Resynthesis of the original signal by $S$ also introduces direct manipulation of the derived synth parameters as an expressive control.

## Approach

A desirable feature of this system is that it be able to operate in real time with minimal latency, i.e. in a live improvisation setting. Therefore we divide the audio stream into short segments and map each segment to a parameter vector. We take a data driven approach to the mapping. Each short windowed segment of audio $a_i$ is converted to a feature vector $b_i = b(a_i)$. Given a training set of feature vectors $B = {b_1, \ldots b_m}$ matched with their corresponding points in the output space $Y = {y_1, \ldots y_m}$, we interpolate $y_i$ from $Y$ by the distances of $b_i$ from elements of $B$. Specifically, we use an average of the $K_{n}$ nearest neighbors weighted distance according to a gaussian:
$$ y_i = \frac{1}{Z}\sum\limits_{j=1}^{K_{n}}exp(-\frac{\lvert b_i-b_j \rvert^2}{2\sigma^2})y_j $$
where $Z$ is a normalizing factor such that the weights sum to 1 and ${b_1, \ldots b_{K_n}}$ are the $K_n$ nearest weights.

## Basis Function

A good basis function should transform audio to some perceptual space. That is, the euclidean distance between examples in feature space should be small iff the corresponding signals are perceived as similar by a human. [Humans decompose audio signals into frequency components, and are deaf to the phase of those components][5]. Therefore, we consider spectral representations of the input signal. We use a [Hann window][4] to reduce artifacts from spectral processing, and an overlap factor of $N_{overlap}$. An easy starting place is a power spectrum of the windowed signal given by the magnitude of the fast fourier transform. [Laroche][6] discusses the meaning of the power spectrum and describes the technique of extracting fine frequency information from the phase spectrum. [Pachet and Aucouturier][3] describe a common basis used for finding timbral similarity between signals, the Mel-frequency cepstrum coefficients. The MFCCs encode a smoothed spectral envelope, which is good for identifying timbre but destroys pitch information.

In the results below, we use the power spectrum weighted to de-emphasize very low and high frequencies.

## Envelope Following

Additionally, we normalize feature vectors to sum to 1, and reapply the scaling factor to the resynthesized signal.

## Synthesis Process

For a synthesis process, we use frequency modulation as described by [Chowning][1]. The process consists of $N_{oscs}$ FM oscillator pairs, each parameterized by a carrier frequency, ratio of modulating frequency to carrier frequency, index of modulation, and weight $c, r, i, w$. The sample value of the process at time $t$ is then: 
$$S(t) = \sum\limits_{i=1}^{N_{oscs}} w \sin( c t + \frac{i}{r}\sin( c r t ))$$
Resynthesized audio is windowed and summed analogously to the analysis windowing.

## Data Generation

We generate data for our model by randomly sampling the output space according to some ad hoc distributions over the space of synth parameters: $y_i \sim P$. Then the corresponding point in feature space is given by $b_i = b(S(y_i))$.

## Hyperparameters

The hyperparameters present can be summarized as following:
- $N_{oscs}$ - number of FM oscillator pairs
- $L_{window}$ - length of analysis window in samples
- $N_{overlap}$ - number of overlapping analysis windows
- $K_{n}$ - number of nearest neighbors to interpolate
- $\sigma$ - width of gaussian for interpolation
- $N_{examples}$ - the size of the generated training set

## Results

The following are a test audio signal and resynthesized versions for a few different hyperparameter values. Images are colored on a log scale. Rows correspond to time windows and columns correspond to features. The sampling rate is 44,100Hz.

<iframe width="676" height="648" scrolling="no" frameborder="no" src="https://w.soundcloud.com/player/?url=https%3A//api.soundcloud.com/tracks/174170819%3Fsecret_token%3Ds-2StCB&amp;auto_play=false&amp;hide_related=false&amp;show_comments=true&amp;show_user=true&amp;show_reposts=false&amp;visual=true"></iframe>

$N_{oscs}=1, L_{window}=2048, N_{overlap}=2, K_{n}=1, N_{examples}=1000$

<iframe width="676" height="648" scrolling="no" frameborder="no" src="https://w.soundcloud.com/player/?url=https%3A//api.soundcloud.com/tracks/174170818&amp;auto_play=false&amp;hide_related=true&amp;show_comments=true&amp;show_user=true&amp;show_reposts=false&amp;visual=true"></iframe>

$N_{oscs}=3, L_{window}=2048, N_{overlap}=2, K_{n}=1, N_{examples}=1000$

<iframe width="676" height="648" scrolling="no" frameborder="no" src="https://w.soundcloud.com/player/?url=https%3A//api.soundcloud.com/tracks/174170812&amp;auto_play=false&amp;hide_related=true&amp;show_comments=true&amp;show_user=true&amp;show_reposts=false&amp;visual=true"></iframe>

$N_{oscs}=3, L_{window}=2048, N_{overlap}=2, K_{n}=1, N_{examples}=10000$

<iframe width="676" height="648" scrolling="no" frameborder="no" src="https://w.soundcloud.com/player/?url=https%3A//api.soundcloud.com/tracks/174170811&amp;auto_play=false&amp;hide_related=true&amp;show_comments=true&amp;show_user=true&amp;show_reposts=false&amp;visual=true"></iframe>

$N_{oscs}=3, L_{window}=2048, N_{overlap}=2, K_{n}=1000, N_{examples}=10000, \sigma = .1$

<iframe width="676" height="648" scrolling="no" frameborder="no" src="https://w.soundcloud.com/player/?url=https%3A//api.soundcloud.com/tracks/174170810&amp;auto_play=false&amp;hide_related=true&amp;show_comments=true&amp;show_user=true&amp;show_reposts=false&amp;visual=true"></iframe>
</div>

The first two resynthesized signals differ by number of oscillators. The second produces more plausible looking spectra, but doesn't really sound any more accurate. Increasing the size of the training set in the third example to compensate the increased dimensionality seems to help. The last example turns on interpolation, which has disastrous results. This is beceause interpolating carrier frequencies and ratios expressed in Hz doesn't make much sense; multiple frequency peaks end up fighting each other, and the result is less accurate than either one. A better representation of synthesis parameters, or perhaps a different synthesis process altogether will be needed to move beyond nearest-neighbor sampling of the training set.


## Future Work

The milestone goal was met. Future work may include:
- Sampling adaptive to a corpus of representative audio to improve quality of training set
- More sophisticated synthesis model, perhaps following [Schottstaedt][2]
- Review the literature on non parametric models
- Continue to develop basis function to represent both spectral envelope and pitch
- Principled tests on effects of hyperparameters

[1]: http://www.jstor.org/stable/23320142 "John M. Chowning. 1977. The Synthesis of Complex Audio Spectra by Means of Frequency Modulation. Computer Music Journal, Vol. 1, No. 2 (April, 1977), pp. 46-54"

[2]: http://www.jstor.org/stable/40731300 "Bill Schottstaedt. The Simulation of Natural Instrument Tones using Frequency Modulation with a Complex Modulating Wave. Computer Music Journal, Vol. 1, No. 4 (November 1977), pp. 46-50"

[3]: http://coltrane.music.mcgill.ca/Andrew/timbre/timbre-similarity_aucoutrier.pdf "Francois Pachet and Jean-Julien Aucouturier. Improving timbre similarity: How high is the sky? Journal of negative results in speech and audio sciences 1.1 (2004): 1-13"

[4]: http://en.wikipedia.org/wiki/Window_function "Wikipedia: Window function"

[5]: http://en.wikipedia.org/wiki/Psychoacoustics "Wikipedia: Psychoacoustics"

[6]: http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=759041 "J. Laroche. Improved phase vocoder time-scale modification of audio. IEEE Transactions on Speech and Audio Processing. Vol. 7, Issue 3 (August 2002), pp 323-332"

