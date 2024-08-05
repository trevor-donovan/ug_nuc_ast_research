Purpose: 
-Integrate 2(E1/2 - E0)n dz to find a nuclear surface tension, and write the answer to console

Properties of the code:
-It uses Simpson's rule (and compares to a Riemann sum) on an even-spaced sampled (not analytic) function. Since each subinterval needs 3 points, there are about half as many subintervals as data points 
-There are some vestigial comments and components from the debugging process which might still come in handy.
-All units should be in MeV and fm.

Included are sample datasets from an SHF code written by Yonsei University's Yeunhwan Lim.