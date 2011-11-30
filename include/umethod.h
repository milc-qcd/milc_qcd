#ifndef _UMETHOD_H
#define _UMETHOD_H

// Unitarization algorithm
// Choices are:
//   UNITARIZE_NONE
//   UNITARIZE_APE -- derivative is not ready
//   UNITARIZE_ROOT -- Y=V*(V^+*V)^-1/2, ()^-1/2 is calculated by inversion,
//                     derivative is calculated with finite difference
//   UNITARIZE_RATIONAL -- Y=V*(V^+*V)^-1/2, ()^-1/2 and its derivative
//                         are calculated with rational approximation
//   UNITARIZE_HISQ -- to be decided
//   UNITARIZE_ANALYTIC -- Y=V*(V^+*V)^-1/2, ()^-1/2 and its derivative
//                         are calculated analytically
//   UNITARIZE_STOUT -- stout smearing from Morningstar, Peardon paper
//   UNITARIZE_ITERATIVE -- iterate V <- 1/2 (V + V^\dagger^-1)
#define UNITARIZE_NONE 0
#define UNITARIZE_APE 1
#define UNITARIZE_ROOT 2
#define UNITARIZE_RATIONAL 3
#define UNITARIZE_HISQ 4
#define UNITARIZE_ANALYTIC 5
#define UNITARIZE_STOUT 6
#define UNITARIZE_ITERATIVE 7

// Unitarization group
// Choices are:
//  UNITARIZE_U3
//  UNITARIZE_SU3
#define UNITARIZE_U3 0
#define UNITARIZE_SU3 1

#endif /* _UMETHOD_H */

