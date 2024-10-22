.SUFFIXES: .test-param .sample-in

.test-param.sample-in::
	build_input.$*.sh $*.test-param > $*.sample-in

TESTIN = \
	ks_spectrum_asqtad.fpi.1.sample-in \
	ks_spectrum_asqtad.nd.1.sample-in \
	ks_spectrum_asqtad.nl.1.sample-in \
	ks_spectrum_asqtad.nlpi2.1.sample-in \
	ks_spectrum_asqtad.spectrum2.1.sample-in \
	ks_spectrum_asqtad.periodic.1.sample-in \
	ks_spectrum_hisq.fpi.1.sample-in \
	ks_spectrum_hisq.nd.1.sample-in \
	ks_spectrum_hisq.nl.1.sample-in \
	ks_spectrum_hisq.nlpi2.1.sample-in \
	ks_spectrum_hisq.spectrum2.1.sample-in \
	ks_spectrum_hisq.periodic.1.sample-in \
	ks_spectrum_asqtad.fpi.2.sample-in \
	ks_spectrum_asqtad.nd.2.sample-in \
	ks_spectrum_asqtad.nl.2.sample-in \
	ks_spectrum_asqtad.nlpi2.2.sample-in \
	ks_spectrum_asqtad.spectrum2.2.sample-in \
	ks_spectrum_asqtad.periodic.2.sample-in \
	ks_spectrum_hisq.fpi.2.sample-in \
	ks_spectrum_hisq.nd.2.sample-in \
	ks_spectrum_hisq.nl.2.sample-in \
	ks_spectrum_hisq.nlpi2.2.sample-in \
	ks_spectrum_hisq.spectrum2.2.sample-in \
	ks_spectrum_hisq.periodic.2.sample-in
	
sample-in: ${TESTIN}


include ../../Make_test_template
