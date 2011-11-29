include ../../Make_test_template

KSFILE1 = ksvector.ans.MdaggerMinv.test.l8
KSFILE2 = ../../binary_samples/ksvector.ans.MdaggerMinv.sample.l8
KSFILE3 = ksvector.ans.Minv.test.l8
KSFILE4 = ../../binary_samples/ksvector.ans.Minv.sample.l8
FFFILE1 = colormatrix.mom.test.l8
FFFILE2 = ../../binary_samples/colormatrix.mom.sample.l8

KSFILES = ${KSFILE1} ${KSFILE2} ${KSFILE3} ${KSFILE4}
FFFILES  = ${FFFILE1} ${FFFILE2}

DIFF_KSVECTOR = ../../file_utilities/diff_ksvector

${DIFF_KSVECTOR}:
	cd ../../file_utilities; make diff_ksvector

DIFF_COLORMATRIX = ../../file_utilities/diff_colormatrix

${DIFF_COLORMATRIX}:
	cd ../../file_utilities; make diff_colormatrix;

ksvector_compare.2: ${KSFILES} ${DIFF_KSVECTOR}
	${LAUNCH} ${DIFF_KSVECTOR} ${KSFILE1} ${KSFILE2} > $@
	${LAUNCH} ${DIFF_KSVECTOR} ${KSFILE3} ${KSFILE4} >> $@
	-/bin/rm ${KSFILE1} ${KSFILE3}

colormatrix_compare.2:  ${FFFILES}  ${DIFF_COLORMATRIX}
	${LAUNCH} ${DIFF_COLORMATRIX}  ${FFFILE1}  ${FFFILE2} > $@
	-/bin/rm ${FFFILE1}

compare_clean:
	-/bin/rm ksvector_compare.2 colormatrix_compare2

test-clean:
	compare_clean



