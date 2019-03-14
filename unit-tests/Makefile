# Author : Firmin Martin
# Last update : 2018 Apr. 1
# Description : Makefile to launch unit tests.

CC=gcc
SRC=../src
INC=-I$(SRC)
# use cmocka for unit testing
LB=-lcmocka -lmpi -lm
# --coverage for lcov
CFLAGS=-g --coverage

BIN= TestUtils TestPerfectHash TestList TestLBList
# required lcov
COV=lcov --capture --directory . --output-file coverage.info -q
GEN_HTML=genhtml coverage.info --output-directory

all: post

pre: 
	$(MAKE) clean

post: main
	$(MAKE) clean

main: pre
	 $(MAKE) --no-print-directory $(BIN)

TestLlist: TestLlist.c $(SRC)/llist.c $(SRC)/reads.c $(SRC)/log.c $(SRC)/mpiMD_utils.c
	$(CC) $^ $(LB) -o $@ $(CFLAGS) $(INC) -D_GNU_SOURCE

TestPerfectHash: TestPerfectHash.c $(SRC)/perfectHash.c $(SRC)/log.c $(SRC)/reads.c $(SRC)/mpiMD_utils.c
	$(CC) $^ $(LB) -lcrypto  -lssl -o $@ $(CFLAGS) $(INC)

TestList: TestList.c $(SRC)/readList.c  $(SRC)/reads.c
	$(CC) $^ $(LB) -o $@ $(CFLAGS) $(INC)

TestLBList: TestLBList.c $(SRC)/readList.c $(SRC)/reads.c $(SRC)/mpiMD_utils.c $(SRC)/log.c
	$(CC) $^ $(LB) -o $@ $(CFLAGS) $(INC)

TestUtils: TestUtils.c $(SRC)/mpiMD_utils.c $(SRC)/log.c
	$(CC) $^ $(LB) -o $@ $(CFLAGS) $(INC)

TestClusterFunctions : TestClusterFunctions.c $(SRC)/log.c $(SRC)/mark_duplicates.c $(SRC)/reads.c $(SRC)/mpiMD_utils.c $(SRC)/readList.c $(SRC)/perfectHash.c $(SRC)/createLBList.c
	$(CC) $^ $(LB) -lcrypto -o $@ $(CFLAGS) $(INC)

.PHONY=clean cov

clean: 
	-rm $(BIN) *.gcda *.gcno


