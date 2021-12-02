SUBDIRS := $(wildcard */.)
SUBDIRS_CLEAN := $(addsuffix .clean, $(SUBDIRS))

all: $(SUBDIRS)

$(SUBDIRS):
	module load openmpi-4.1.1+gnu-9.3.0
	$(MAKE) -C $@

clean: $(SUBDIRS_CLEAN)

$(SUBDIRS_CLEAN):
	$(MAKE) -C $(basename $@) clean


.PHONY: all $(SUBDIRS) clean $(SUBDIRS_CLEAN)
