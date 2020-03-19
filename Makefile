# Pipeline to build Penalized time-varying COX model

current: target
-include target.mk

######################################################################

ms = makestuff
Sources += $(wildcard *.R *.Rmd)
Sources += Makefile rmd.mk

######################################################################


## Time invariant coxph

### Soft-thresholding operater
proxupdate.Rout: proxupdate.R

### Negative log likelihood
nloglik.Rout: nloglik.R

### Proximal gradient
gradient.Rout: gradient.R

### Other helper functions
helperfuns.Rout: helperfuns.R

### Penalized cox frame
pencox.Rout: pencox.R

### Writeup
pencox.pdf: pencox.rmd

### Example
examples.Rout: examples.R


######################################################################

clean: 
	rm -f *Rout.*  *.Rout .*.RData .*.Rout.* .*.wrapR.* .*.Rlog *.RData *.wrapR.* *.Rlog

######################################################################

### Makestuff

Ignore += makestuff
msrepo = https://github.com/dushoff
Makefile: makestuff/Makefile
makestuff/Makefile:
	git clone $(msrepo)/makestuff
	ls $@

include rmd.mk
-include makestuff/os.mk
-include makestuff/visual.mk
-include makestuff/projdir.mk
-include makestuff/texdeps.mk
-include makestuff/pandoc.mk
-include makestuff/stepR.mk
-include makestuff/git.mk

