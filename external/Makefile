.PHONY: futils distclean make_link lib
.ONESHELL:

FUTILS_TO_USE=futils-gene

# order does matter, first make_link then look for the library
futils: make_link
	@make lib

make_link: $(FUTILS_TO_USE)
	@if [ ! -d $(MACHINE) ] ; then mkdir $(MACHINE); fi
	@if [ -L $(MACHINE)/futils ] ; then rm $(MACHINE)/futils; fi
	@if [ -d futils ] ; then mv futils futils_old; fi
	@if [ ! -d $(MACHINE)/$(FUTILS_TO_USE) ] ; then \
		if [ -d $(FUTILS_TO_USE) ]; then \
			cp -R $(FUTILS_TO_USE) $(MACHINE); \
			ln -s -f $(FUTILS_TO_USE) $(MACHINE)/futils; \
		fi; \
	else \
		ln -s -f $(FUTILS_TO_USE) $(MACHINE)/futils; \
	fi


futils-gene.tar.gz:
	@(wget http://documents.epfl.ch/users/t/tt/ttran/www/$@) || \
	(echo; echo \*\*\* cannot connect to http://documents.epfl.ch/users/t/tt/ttran/www/$@ \*\*\*;\
	echo \*\*\* please try to copy $@ by hand to $(EXTDIR) \*\*\*; echo;)

futils-gene.tar: futils-gene.tar.gz
	gunzip -f futils-gene.tar.gz

futils-gene: futils-gene.tar
	tar -x -f futils-gene.tar
	cp -f Makefile.$(FUTILS_TO_USE) $(FUTILS_TO_USE)/src/Makefile

futils-2.0.tar.gz:
	(wget http://documents.epfl.ch/users/t/tt/ttran/www/$@) || \
	(echo; echo \*\*\* cannot connect to http://documents.epfl.ch/users/t/tt/ttran/www/$@ \*\*\*;\
	echo \*\*\* please try to copy $@ by hand to $(EXTDIR) \*\*\*; echo;)

futils-2.0.tar: futils-2.0.tar.gz
	gunzip -f futils-2.0.tar.gz

futils-2.0: futils-2.0.tar
	tar --extract --strip-components=1 -f futils-2.0.tar

lib: $(FUTILSDIR)/libfutils.a

$(FUTILSDIR)/libfutils.a: $(FUTILSDIR)/futils.f90
	@if [ -s $(FUTILSDIR)/Makefile.$(MACHINE) ] ;\
	then FUTILSMAKEFILE=Makefile.$(MACHINE) ;\
	else FUTILSMAKEFILE=Makefile ;\
	fi ;\
	echo "futilsmakefile = $$FUTILSMAKEFILE";\
	make "F90=$(MPFC)" "FFLAGS=$(FFLAGS)" \
	"LD_FLAGS=$(LDFLAGS)" "C_FLAGS=$(CFLAGS)" \
	"HDF5LIBS=$(HDF5_LIBS)" "HDF5=$(HDF5PATH)" \
	"AR=$(ARCHIVE)" -C $(FUTILSDIR) \
	-f $$FUTILSMAKEFILE lib

distclean:
	@rm -f futils-gene.tar* futils-2.0.tar*
	@rm -rf futils-gene $(MACHINE)

#	@make -C $(FUTILSDIR) clean