# speechfiltering
A MATLAB script to filter speech from 2 channel, uncompressed, .WAV files.

Filtering is done in 3 types (direct convolution, fft convolution, overlap add convolution)

Their speeds are compared, and are relatively fast. The FFT Convolution performs at a rate of 11.949 GSPS!!!

Requirements: MATLAB r2018 (any other version should work, no special libraries needed), and a .wav File longer than 1 min 8 seconds.

Please see attached powerpoint presentation (PDF).
