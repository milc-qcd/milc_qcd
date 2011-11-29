include ../../Make_test_template

W_FILELIST = w_source.1.test-out ks_source_naive.1.test-out
KS_FILELIST = ks_source.1.test-out

# To test the source files the code generates, we don't compare binary
# files, but we do compare checksums

src.1.test-out: ${W_FILELIST} ${KS_FILELIST}
	echo "BEGIN" > $@; \
	../get_checksums_ks.sh ${KS_FILELIST} >> $@;\
	../get_checksums_w.sh  ${W_FILELIST}  >> $@; \
	echo "END" >> $@


