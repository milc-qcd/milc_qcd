include ../../Make_test_template

KSFILE1 = ksprop_scidac1.test.l4448.ll
KSFILE2 = ksprop_scidac2.test.l4448.ll
KSFILE3 = ksprop_scidac3.test.l4448.ll

WFILE1 = wprop_fnal1.test.l4448.ll
WFILE2 = wprop_fnal2.test.l4448.ll
WFILE3 = wprop_fnal3.test.l4448.ll

KSFILES = ${KSFILE1} ${KSFILE2}
WFILES  = ${WFILE1} ${WFILE2}

../diff_ksprop:
	cd .. ; make diff_ksprop

# Usually we don't look at stdout, but we see stderr from diff_wprop
ksprop.test.compare: ${KSFILES} ../diff_ksprop
	${LAUNCH} ../diff_ksprop ${KSFILE1} ${KSFILE2} > ksprop.test.compare
	${LAUNCH} ../diff_ksprop ${KSFILE1} ${KSFILE2} >> ksprop.test.compare

../diff_wprop:
	cd .. ; make diff_wprop

# Usually we don't look at stdout, but we see stderr from diff_wprop
wprop.test.compare:  ${WFILES} ../diff_wprop
	${LAUNCH} ../diff_wprop  ${WFILE1}  ${WFILE2} > wprop.test.compare
	${LAUNCH} ../diff_wprop  ${WFILE1}  ${WFILE2} >> wprop.test.compare

compare_clean:
	-/bin/rm ksprop.test.compare wprop.test.compare ${KSFILES} ${WFILES}

test-clean:
	compare_clean
