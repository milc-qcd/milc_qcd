SU3 Wilson valence fermions;  PROPAGATING form factor code 
Heavy--> Heavy and Heavy-->light form factors
MIMD version 5
Using CLOVER fermions with the BI-CG and hopping algorithms
Machine = Scalar processor, with 1 nodes
start: Fri Apr 30 20:51:41 1999

type 0 for no prompts  or 1 for prompts
nx 4
ny 4
nz 4
nt 8
LAYOUT = Timeslices, options = EVENFIRST,
Mallocing 3.0 MBytes per node for lattice


verbose_flag 0
clov_c 1.1
u0 0.836
reload_serial lat.sample.l4448
nkap_spectator 1
max_cg_iterations 100
max_cg_restarts 10
error_for_propagator 1e-05
kappa_spectator 0.12
gaussian
r0 2
reload_serial_prop lprop.sample.l4448
nkap_light_zonked 1
max_cg_iterations 100
max_cg_restarts 10
error_for_propagator 1e-05
kappa_zonked_light 0.12
gaussian
r0 2
reload_serial_prop lprop.sample.l4448
save_serial_prop _ssink
nkap_heavy_zonked 1
max_cg_iterations 100
max_cg_restarts 10
error_for_propagator 1e-05
kappa_zonked_heavy 0.12
hopilu
gaussian
r0 2
save_serial_prop hprop.test.l4448
save_serial_prop _ssink
nkap_sequential 1
kappa_seq 0.12
hopilu
final_time 3
start_of_momentum p
p= 0   0   0
start_of_momentum q
q= 0   0   0
start_of_momentum k
k= 0   0   0
seq_smear_func in.sample.local_smear
save_binary HH3 out.test.HH3_l4448
save_binary HL3 out.test.HL3_l4448
save_ascii HH2_GL out.test.HH2_GL_l4448
save_binary LL2_GG out.test.LL2_GG_l4448
save_binary HL2_GG out.test.HL2_GG_l4448
save_binary HL2_GE out.test.HL2_GE_l4448
save_binary HL2_GL out.test.HL2_GL_l4448
Reading as old-style gauge field configuration.
Old format header parameters are 5.500000 0.050000
Restored binary gauge configuration serially from file lat.sample.l4448
Time to reload gauge configuration = 1.000000e-02
CHECK PLAQ: 1.723750e+00 1.690586e+00
Unitarity checked.  Max deviation 2.38e-07
Time to check unitarity = 1.000000e-02
load_scalar_smear: I have read a smearing function from the file in.sample.local_smear
Allocating 0.0 MBytes per node for HL3 correlator
Allocating 0.0 MBytes per node for HH3 correlator
Allocating 0.0 MBytes per node for LL2 correlator
Allocating 0.0 MBytes per node for HL2 correlator
Allocating 0.0 MBytes per node for HL2 correlator

STARTING SMEARED COLOR 0
REUSING previous zonked light tmp file lprop.sample.l4448_ssink
Time to reload prop spin 0 color 0 0.000000e+00
Time to reload prop spin 1 color 0 1.000000e-02
Time to reload prop spin 2 color 0 0.000000e+00
Time to reload prop spin 3 color 0 0.000000e+00
Starting heavy_zonked inversions, kappa = 0.120000
HOPILU: time = 1.23e+00 size_r= 9.80e-06 iters= 32 MF = 24.9
 OK converged size_r= 9.8e-06 iters= 32
Starting heavy_zonked inversions, kappa = 0.120000
HOPILU: time = 1.28e+00 size_r= 8.11e-06 iters= 32 MF = 24.0
 OK converged size_r= 8.1e-06 iters= 32
Starting heavy_zonked inversions, kappa = 0.120000
HOPILU: time = 1.25e+00 size_r= 9.00e-06 iters= 32 MF = 24.5
 OK converged size_r= 9e-06 iters= 32
Starting heavy_zonked inversions, kappa = 0.120000
HOPILU: time = 1.28e+00 size_r= 8.20e-06 iters= 33 MF = 24.7
 OK converged size_r= 8.2e-06 iters= 33
Time to save prop spin 0 color 0 = 1.000000e-02
Time to save prop spin 1 color 0 = 1.000000e-02
Time to save prop spin 2 color 0 = 0.000000e+00
Time to save prop spin 3 color 0 = 1.000000e-02
Time to save prop spin 0 color 0 = 0.000000e+00
Time to save prop spin 1 color 0 = 1.000000e-02
Time to save prop spin 2 color 0 = 1.000000e-02
Time to save prop spin 3 color 0 = 1.000000e-02
Computing LL2_GG correlator
meson_cont_mom_lean2: time 0.2 sec 19.5 MF
Computing HL2_GG correlator
meson_cont_mom_lean2: time 0.2 sec 18.7 MF

STARTING SMEARED COLOR 1
REUSING previous zonked light tmp file lprop.sample.l4448_ssink
Time to reload prop spin 0 color 1 0.000000e+00
Time to reload prop spin 1 color 1 1.000000e-02
Time to reload prop spin 2 color 1 1.000000e-02
Time to reload prop spin 3 color 1 0.000000e+00
Starting heavy_zonked inversions, kappa = 0.120000
HOPILU: time = 1.31e+00 size_r= 8.58e-06 iters= 33 MF = 24.1
 OK converged size_r= 8.6e-06 iters= 33
Starting heavy_zonked inversions, kappa = 0.120000
HOPILU: time = 1.33e+00 size_r= 7.75e-06 iters= 34 MF = 24.5
 OK converged size_r= 7.8e-06 iters= 34
Starting heavy_zonked inversions, kappa = 0.120000
HOPILU: time = 1.30e+00 size_r= 7.45e-06 iters= 34 MF = 25.1
 OK converged size_r= 7.4e-06 iters= 34
Starting heavy_zonked inversions, kappa = 0.120000
HOPILU: time = 1.26e+00 size_r= 7.58e-06 iters= 33 MF = 25.1
 OK converged size_r= 7.6e-06 iters= 33
Time to save prop spin 0 color 1 = 1.000000e-02
Time to save prop spin 1 color 1 = 0.000000e+00
Time to save prop spin 2 color 1 = 1.000000e-02
Time to save prop spin 3 color 1 = 1.000000e-02
Time to save prop spin 0 color 1 = 1.000000e-02
Time to save prop spin 1 color 1 = 1.000000e-02
Time to save prop spin 2 color 1 = 0.000000e+00
Time to save prop spin 3 color 1 = 1.000000e-02
Computing LL2_GG correlator
meson_cont_mom_lean2: time 0.2 sec 17.9 MF
Computing HL2_GG correlator
meson_cont_mom_lean2: time 0.2 sec 18.7 MF

STARTING SMEARED COLOR 2
REUSING previous zonked light tmp file lprop.sample.l4448_ssink
Time to reload prop spin 0 color 2 0.000000e+00
Time to reload prop spin 1 color 2 0.000000e+00
Time to reload prop spin 2 color 2 1.000000e-02
Time to reload prop spin 3 color 2 1.000000e-02
Starting heavy_zonked inversions, kappa = 0.120000
HOPILU: time = 1.31e+00 size_r= 7.49e-06 iters= 33 MF = 24.1
 OK converged size_r= 7.5e-06 iters= 33
Starting heavy_zonked inversions, kappa = 0.120000
HOPILU: time = 1.23e+00 size_r= 8.39e-06 iters= 32 MF = 24.9
 OK converged size_r= 8.4e-06 iters= 32
Starting heavy_zonked inversions, kappa = 0.120000
HOPILU: time = 1.26e+00 size_r= 9.04e-06 iters= 33 MF = 25.1
 OK converged size_r= 9e-06 iters= 33
Starting heavy_zonked inversions, kappa = 0.120000
HOPILU: time = 1.30e+00 size_r= 9.20e-06 iters= 33 MF = 24.3
 OK converged size_r= 9.2e-06 iters= 33
Time to save prop spin 0 color 2 = 0.000000e+00
Time to save prop spin 1 color 2 = 1.000000e-02
Time to save prop spin 2 color 2 = 0.000000e+00
Time to save prop spin 3 color 2 = 1.000000e-02
Time to save prop spin 0 color 2 = 0.000000e+00
Time to save prop spin 1 color 2 = 1.000000e-02
Time to save prop spin 2 color 2 = 1.000000e-02
Time to save prop spin 3 color 2 = 1.000000e-02
Computing LL2_GG correlator
meson_cont_mom_lean2: time 0.2 sec 19.5 MF
Computing HL2_GG correlator
meson_cont_mom_lean2: time 0.2 sec 17.9 MF
I have written the two point functions to the file out.test.LL2_GG_l4448
I have written the two point functions to the file out.test.HL2_GG_l4448
Wrote prop file hprop.test.l4448 time stamp Fri Apr 30 20:51:41 1999
Wrote prop file hprop.test.l4448_ssink time stamp Fri Apr 30 20:51:41 1999

STARTING LOCAL COLOR 0
Reading as 1996-style propagator field configuration.
Time to reload prop spin 0 color 0 0.000000e+00
BiCGILU: time = 1.14e+00 size_r= 8.82e-06 iters= 13 MF = 25.5
 OK converged size_r= 8.8e-06 iters= 13
Reading as 1996-style propagator field configuration.
Time to reload prop spin 1 color 0 1.000000e-02
BiCGILU: time = 1.15e+00 size_r= 4.42e-06 iters= 13 MF = 25.3
 OK converged size_r= 4.4e-06 iters= 13
Reading as 1996-style propagator field configuration.
Time to reload prop spin 2 color 0 1.000000e-02
BiCGILU: time = 1.14e+00 size_r= 9.76e-06 iters= 13 MF = 25.5
 OK converged size_r= 9.8e-06 iters= 13
Reading as 1996-style propagator field configuration.
Time to reload prop spin 3 color 0 1.000000e-02
BiCGILU: time = 1.28e+00 size_r= 5.51e-06 iters= 14 MF = 24.5
 OK converged size_r= 5.5e-06 iters= 14
Time to reload prop spin 0 color 0 0.000000e+00
Time to reload prop spin 1 color 0 1.000000e-02
Time to reload prop spin 2 color 0 1.000000e-02
Time to reload prop spin 3 color 0 1.000000e-02
Computing HL2_GL correlator
meson_cont_mom_lean2: time 0.2 sec 17.9 MF
meson_cont_mom_lean2: time 0.2 sec 17.9 MF
meson_cont_mom_lean2: time 0.2 sec 18.7 MF

Starting momentum insertion 0
Sequential smearing for momentum insertion 0
Computing HL2_GE correlator
meson_cont_mom_lean2: time 0.2 sec 19.5 MF
Computing form factors for k_sequential = 0.12 k_spectator = 0.12 p = 0,0,0 
Sequential source inversion Kappa = 0.120000 
HOPILU: time = 1.52e+00 size_r= 7.54e-06 iters= 40 MF = 25.2
 OK converged size_r= 7.5e-06 iters= 40
Sequential source inversion Kappa = 0.120000 
HOPILU: time = 1.46e+00 size_r= 9.82e-06 iters= 38 MF = 24.9
 OK converged size_r= 9.8e-06 iters= 38
Sequential source inversion Kappa = 0.120000 
HOPILU: time = 1.56e+00 size_r= 9.84e-06 iters= 39 MF = 24.0
 OK converged size_r= 9.8e-06 iters= 39
Sequential source inversion Kappa = 0.120000 
HOPILU: time = 1.46e+00 size_r= 9.67e-06 iters= 38 MF = 24.9
 OK converged size_r= 9.7e-06 iters= 38
meson_cont_mom_lean2: time 0.2 sec 18.7 MF
meson_cont_mom_lean2: time 0.2 sec 18.7 MF
meson_cont_mom_lean2: time 0.2 sec 18.7 MF
meson_cont_mom_lean2: time 0.2 sec 19.5 MF
meson_cont_mom_lean2: time 0.2 sec 18.7 MF
meson_cont_mom_lean2: time 0.2 sec 19.5 MF

STARTING LOCAL COLOR 1
Reading as 1996-style propagator field configuration.
Time to reload prop spin 0 color 1 1.000000e-02
BiCGILU: time = 1.25e+00 size_r= 7.95e-06 iters= 14 MF = 25.1
 OK converged size_r= 7.9e-06 iters= 14
Reading as 1996-style propagator field configuration.
Time to reload prop spin 1 color 1 1.000000e-02
BiCGILU: time = 1.46e+00 size_r= 1.59e-06 iters= 16 MF = 24.5
 OK converged size_r= 1.6e-06 iters= 16
Reading as 1996-style propagator field configuration.
Time to reload prop spin 2 color 1 0.000000e+00
BiCGILU: time = 1.23e+00 size_r= 4.85e-06 iters= 14 MF = 25.5
 OK converged size_r= 4.9e-06 iters= 14
Reading as 1996-style propagator field configuration.
Time to reload prop spin 3 color 1 1.000000e-02
BiCGILU: time = 1.23e+00 size_r= 4.83e-06 iters= 14 MF = 25.5
 OK converged size_r= 4.8e-06 iters= 14
Time to reload prop spin 0 color 1 0.000000e+00
Time to reload prop spin 1 color 1 1.000000e-02
Time to reload prop spin 2 color 1 1.000000e-02
Time to reload prop spin 3 color 1 1.000000e-02
Computing HL2_GL correlator
meson_cont_mom_lean2: time 0.2 sec 19.5 MF
meson_cont_mom_lean2: time 0.2 sec 17.9 MF
meson_cont_mom_lean2: time 0.2 sec 19.5 MF

Starting momentum insertion 0
Sequential smearing for momentum insertion 0
Computing HL2_GE correlator
meson_cont_mom_lean2: time 0.2 sec 18.7 MF
Computing form factors for k_sequential = 0.12 k_spectator = 0.12 p = 0,0,0 
Sequential source inversion Kappa = 0.120000 
HOPILU: time = 1.55e+00 size_r= 9.37e-06 iters= 39 MF = 24.1
 OK converged size_r= 9.4e-06 iters= 39
Sequential source inversion Kappa = 0.120000 
HOPILU: time = 1.50e+00 size_r= 8.45e-06 iters= 39 MF = 24.9
 OK converged size_r= 8.5e-06 iters= 39
Sequential source inversion Kappa = 0.120000 
HOPILU: time = 1.49e+00 size_r= 9.41e-06 iters= 39 MF = 25.1
 OK converged size_r= 9.4e-06 iters= 39
Sequential source inversion Kappa = 0.120000 
HOPILU: time = 1.48e+00 size_r= 8.67e-06 iters= 38 MF = 24.6
 OK converged size_r= 8.7e-06 iters= 38
meson_cont_mom_lean2: time 0.2 sec 18.7 MF
meson_cont_mom_lean2: time 0.2 sec 17.9 MF
meson_cont_mom_lean2: time 0.2 sec 18.7 MF
meson_cont_mom_lean2: time 0.2 sec 17.9 MF
meson_cont_mom_lean2: time 0.2 sec 17.9 MF
meson_cont_mom_lean2: time 0.2 sec 19.5 MF

STARTING LOCAL COLOR 2
Reading as 1996-style propagator field configuration.
Time to reload prop spin 0 color 2 1.000000e-02
BiCGILU: time = 1.15e+00 size_r= 5.36e-06 iters= 13 MF = 25.3
 OK converged size_r= 5.4e-06 iters= 13
Reading as 1996-style propagator field configuration.
Time to reload prop spin 1 color 2 1.000000e-02
BiCGILU: time = 1.15e+00 size_r= 9.82e-06 iters= 13 MF = 25.3
 OK converged size_r= 9.8e-06 iters= 13
Reading as 1996-style propagator field configuration.
Time to reload prop spin 2 color 2 0.000000e+00
BiCGILU: time = 1.25e+00 size_r= 2.74e-06 iters= 14 MF = 25.1
 OK converged size_r= 2.7e-06 iters= 14
Reading as 1996-style propagator field configuration.
Time to reload prop spin 3 color 2 0.000000e+00
BiCGILU: time = 1.21e+00 size_r= 9.21e-06 iters= 13 MF = 24.0
 OK converged size_r= 9.2e-06 iters= 13
Time to reload prop spin 0 color 2 1.000000e-02
Time to reload prop spin 1 color 2 0.000000e+00
Time to reload prop spin 2 color 2 1.000000e-02
Time to reload prop spin 3 color 2 1.000000e-02
Computing HL2_GL correlator
meson_cont_mom_lean2: time 0.2 sec 18.7 MF
meson_cont_mom_lean2: time 0.2 sec 18.7 MF
meson_cont_mom_lean2: time 0.2 sec 19.5 MF

Starting momentum insertion 0
Sequential smearing for momentum insertion 0
Computing HL2_GE correlator
meson_cont_mom_lean2: time 0.2 sec 18.7 MF
Computing form factors for k_sequential = 0.12 k_spectator = 0.12 p = 0,0,0 
Sequential source inversion Kappa = 0.120000 
HOPILU: time = 1.49e+00 size_r= 8.91e-06 iters= 39 MF = 25.1
 OK converged size_r= 8.9e-06 iters= 39
Sequential source inversion Kappa = 0.120000 
HOPILU: time = 1.42e+00 size_r= 9.43e-06 iters= 37 MF = 25.0
 OK converged size_r= 9.4e-06 iters= 37
Sequential source inversion Kappa = 0.120000 
HOPILU: time = 1.55e+00 size_r= 8.71e-06 iters= 39 MF = 24.1
 OK converged size_r= 8.7e-06 iters= 39
Sequential source inversion Kappa = 0.120000 
HOPILU: time = 1.49e+00 size_r= 7.95e-06 iters= 39 MF = 25.1
 OK converged size_r= 8e-06 iters= 39
meson_cont_mom_lean2: time 0.2 sec 19.5 MF
meson_cont_mom_lean2: time 0.2 sec 18.7 MF
meson_cont_mom_lean2: time 0.2 sec 18.7 MF
meson_cont_mom_lean2: time 0.2 sec 19.5 MF
meson_cont_mom_lean2: time 0.2 sec 18.7 MF
meson_cont_mom_lean2: time 0.2 sec 18.7 MF
I have written propagating form factors to the file out.test.HL3_l4448
I have written propagating form factors to the file out.test.HH3_l4448
I have written the two point functions to the file out.test.HL2_GE_l4448
I have written the two point functions to the file out.test.HL2_GL_l4448
RUNNING COMPLETED
Total time = 68.86 sec


exit: Fri Apr 30 20:52:51 1999

