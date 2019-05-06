
bin_PROGRAMS = mpiMD
noinst_HEADERS = $(mpiMDHDR) $(mpiSortHDR)

mpiMD_SOURCES = $(mpiMDSRC) $(mpiSortSRC)

mpiMDSRC = mark_duplicates.c mark_duplicates_discordant.c createLBList.c perfectHash.c log.c mpiMD_utils.c reads.c llist.c

mpiMDHDR = mark_duplicates.h mark_duplicates_discordant.h createLBList.h perfectHash.h log.h mpiMD_utils.h common.h reads.h llist.h

mpiSortSRC = mpiSort.c mpiSort_utils.c mergeSort.c parabitonicsort.c \
	parabitonicsort2.c parabitonicsort3.c qksort.c parser.c \
	tokenizer.c preWrite.c write.c  write_utils.c sort_any_dim.c

mpiSortHDR = mpiSort_utils.h mergeSort.h parabitonicsort.h \
	parabitonicsort2.h parabitonicsort3.h qksort.h \
	parser.h tokenizer.h preWrite.h write.h write_utils.h \
	bgzf.c bgzf.h khash.h sort_any_dim.h \
	compat.h

AM_CFLAGS= -std=gnu99 
#AM_LDFLAGS=--coverage
