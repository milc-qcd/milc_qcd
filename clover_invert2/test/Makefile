../tie_open_meson:
	cd .. ; make tie_open_meson

../dump_open_meson:
	cd .. ; make dump_open_meson

tie_open_meson.1.test-out: open_meson.1.test-out ../tie_open_meson
	../tie_open_meson < tie_open_meson.1.sample-in > tie_open_meson.1.test-out

dump_open_meson.1.test-out: open_meson.1.test-out ../dump_open_meson
	../dump_open_meson < dump_open_meson.1.sample-in > dump_open_meson.1.test-out

tie_open_meson.2.test-out: open_meson.2.test-out ../tie_open_meson
	../tie_open_meson < tie_open_meson.2.sample-in > tie_open_meson.2.test-out

dump_open_meson.2.test-out: open_meson.2.test-out ../dump_open_meson
	../dump_open_meson < dump_open_meson.2.sample-in > dump_open_meson.2.test-out

tie_open_meson.1.sample-out: open_meson.1.sample-out ../tie_open_meson
	../tie_open_meson < tie_open_meson.1.sample-in > tie_open_meson.1.sample-out

dump_open_meson.1.sample-out: open_meson.1.sample-out ../dump_open_meson
	../dump_open_meson < dump_open_meson.1.sample-in > dump_open_meson.1.sample-out

tie_open_meson.2.sample-out: open_meson.2.sample-out ../tie_open_meson
	../tie_open_meson < tie_open_meson.2.sample-in > tie_open_meson.2.sample-out

dump_open_meson.2.sample-out: open_meson.2.sample-out ../dump_open_meson
	../dump_open_meson < dump_open_meson.2.sample-in > dump_open_meson.2.sample-out

include ../../Make_test_template

