    prompt 0
    nx 4
    ny 4
    nz 4
    nt 8

   verbose_flag 0
   reload_serial  ../binary_samples/lat.sample.l4448

   nkap_spectator  1
   max_cg_iterations 100
   max_cg_restarts 10
   error_for_propagator .00001
   gaussian
   kappa_spectator  0.120
   r0 1.0
   fresh_wprop

   nkap_light_zonked 1
   max_cg_iterations 100
   max_cg_restarts 10
   error_for_propagator .00001
   gaussian
   kappa_zonked_light 0.120
   r0 1.0
   fresh_wprop

   vel_file in.sample.input_vel
   mom_file  in.sample.input_mom
   heavy_light_out  out.test.form_out_l4448_w
   twopt_out  out.test.two_out_l4448_w
   seq_out  out.test.seq_out_l4448_w


