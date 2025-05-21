# PHAFD (parallel hybrid atom-field dynamics)


From Discrete Fourier Transform, we have that

f(q) = L^3/N^3*DFT[f(x)]
f(x) = 1/L^3 IDFT[f(q)]

so that to normalise we have

f(x) = 1/N^3 IDFT[DFT[f(x)]],

i.e. a factor of 1/N^3 must be included at each step. We choose to
put this factor in right before the final inverse fourier transform.