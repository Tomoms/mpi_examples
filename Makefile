SUBDIRS := $(wildcard */.)
SUBDIRS_CLEAN := $(addsuffix .clean, $(SUBDIRS))

all: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

clean: $(SUBDIRS_CLEAN)

$(SUBDIRS_CLEAN):
	$(MAKE) -C $(basename $@) clean


.PHONY: all $(SUBDIRS) clean $(SUBDIRS_CLEAN)
