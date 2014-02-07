SRCS = ccamcoverage.eps ccameotrace.eps ccam.eps ccamgnatrace.eps ccamimstats.eps ccamiotrace.eps ccamopstats.eps ccamoptrace.eps ccamx0.eps ccamxfinal.eps roma.eps

PDFS	= ${SRCS:.eps=.pdf}

# clear out all suffixes
.SUFFIXES:
# list only those we use 
.SUFFIXES: .eps .pdf

# define a pattern rule for .eps -> .pdf
%.pdf : %.eps
	ps2pdf -dEPSCrop $< $@

all: $(PDFS)
