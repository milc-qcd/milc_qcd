include ../../Make_test_template

.SUFFIXES: .out .test-out

OUTLIST1 = \
	O5_ll.m0.010000.1.out \
	O5_ll.m0.020000.1.out \
	O5_ls.m0.010000.s0.050000.1.out \
	O5_ls.m0.020000.s0.050000.1.out \
	Omu_lHH.m0.010000.1.out \
	Omu_lHH.m0.020000.1.out \
	Omu_ll.m0.010000.1.out \
	Omu_ll.m0.020000.1.out \
	Omu_ls.m0.010000.s0.050000.1.out \
	Omu_ls.m0.020000.s0.050000.1.out \
	Omu_sHH.s0.050000.1.out \
	Omu_ss.s0.050000.1.out

OUTLIST2 = \
	O5_ll.m0.010000.2.out \
	O5_ll.m0.020000.2.out \
	O5_ls.m0.010000.s0.050000.2.out \
	O5_ls.m0.020000.s0.050000.2.out \
	Omu_lHH.m0.010000.2.out \
	Omu_lHH.m0.020000.2.out \
	Omu_ll.m0.010000.2.out \
	Omu_ll.m0.020000.2.out \
	Omu_ls.m0.010000.s0.050000.2.out \
	Omu_ls.m0.020000.s0.050000.2.out \
	Omu_sHH.s0.050000.2.out \
	Omu_ss.s0.050000.2.out

${OUTLIST1}: ks_hl_baryon.1.test-out

${OUTLIST2}: ks_hl_baryon.2.test-out

.out.test-out::
	/bin/mv $*.out $*.test-out
