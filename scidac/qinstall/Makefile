# Library versions and git branch tags

QMPV= quda-1.1.0-consistent
QMPT= 3010fef  

QIOV= quda-1.1.0-consistent
QIOT= a5c3ae580b846130c06dc060ab822f17d6fe2171

QLAV= 1.9.0
QLAT= qla1-9-0

QDPV= 1.11.1
QDPT= qdp1-11-1

#QOPV = 0.21.1
#QOPT = qopqdp0-21-1

QOPV = master
QOPT = master

install: install_qmp install_qio install_qla install_qdp  install_qopqdp

install_qmp: ../install/qmp-$(QMPV)-openmpi

install_qio: ../install/qio-$(QIOV)-openmpi

install_qla: ../install/qla-$(QLAV)

install_qdp: ../install/qdp-$(QDPV)

install_qopqdp: ../install/qopqdp-$(QOPV)

# QMP targets

../src/qmp-$(QMPV):
	@echo Checkout $(QMPT)
	( git clone https://github.com/usqcd-software/qmp.git ../src/qmp-$(QMPV) && \
	  cd ../src/qmp-$(QMPV) && \
          git checkout $(QMPT) && \
          autoreconf -fi )

../install/qmp-$(QMPV)-single: ../src/qmp-$(QMPV)
	@echo Install qmp-$(QMPV)-single
	./qinstall generic-single qmp $(QMPV)

../install/qmp-$(QMPV)-openmpi: ../src/qmp-$(QMPV)
	@echo Install qmp-$(QMPV)-openmpi
	./qinstall generic qmp $(QMPV)

# QIO targets

../src/qio-$(QIOV):
	@echo Checkout $(QIOT)
	( git clone https://github.com/usqcd-software/qio.git ../src/qio-$(QIOV) && \
	  cd ../src/qio-$(QIOV) && \
	  git checkout $(QIOT) && \
	  git clone https://github.com/usqcd-software/c-lime.git other_libs/c-lime && \
	  ( cd other_libs/c-lime && autoreconf -fi ) && autoreconf -fi )

../install/qio-$(QIOV)-single: ../install/qmp-$(QMPV)-single ../src/qio-$(QIOV)
	@echo Install qio-$(QIOV)-single
	./qinstall generic-single qio $(QIOV)

../install/qio-$(QIOV)-openmpi: ../install/qmp-$(QMPV)-openmpi ../src/qio-$(QIOV)
	@echo Install qio-$(QIOV)-openmpi
	./qinstall generic qio $(QIOV)

# QLA targets

../src/qla-$(QLAV):
	@echo Checkout $(QLAT)
	git clone https://github.com/usqcd-software/qla.git --branch $(QLAT) --single-branch ../src/qla-$(QLAV)
	( cd ../src/qla-$(QLAV) && autoreconf -fi )

../install/qla-$(QLAV): ../src/qla-$(QLAV)
	@echo Install qla-$(QLAV)
	./qinstall generic qla $(QLAV)

# QDP targets

../src/qdp-$(QDPV):
	@echo Checkout $(QDPT)
	git clone https://github.com/usqcd-software/qdp.git --branch $(QDPT) --single-branch ../src/qdp-$(QDPV)
	( cd ../src/qdp-$(QDPV) && autoreconf -fi )

../install/qdp-$(QDPV): ../install/qmp-$(QMPV)-single ../install/qla-$(QLAV) ../src/qdp-$(QDPV)
	@echo Install qdp-$(QDPV)
	./qinstall generic-single qdp $(QDPV)

# QOPQDP targets
# TODO: revert to a target that produces a qopqdp.a file; Add target(s) to also build partial-linked dynamic libs: one for single, one linked against MPI?

../src/qopqdp-$(QOPV):
	@echo Checkout $(QOPT)
	git clone https://github.com/usqcd-software/qopqdp.git --branch $(QOPT) --single-branch ../src/qopqdp-$(QOPV)
	( cd ../src/qopqdp-$(QOPV) && autoreconf -fi )

../install/qopqdp-$(QOPV): ../install/qio-$(QIOV)-single ../install/qdp-$(QDPV) ../src/qopqdp-$(QOPV)
	@echo Install qopqdp-$(QOPV)
	./qinstall generic-single qopqdp $(QOPV)

