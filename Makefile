CHILD_PATH = .
SUBPATHS = $(shell find $(CHILD_PATH) -mindepth 1 -maxdepth 1 -name "*" -type d | sort)

.DEFAULT_GOAL : all
all :
	@for i in $(SUBPATHS); do \
	echo "make $@ in $$i..."; \
	(cd $$i; $(MAKE) $@); done

.PHONY : clean
clean :
	@for i in $(SUBPATHS); do \
	echo "make $@ in $$i..."; \
	(cd $$i; $(MAKE) $@); done
