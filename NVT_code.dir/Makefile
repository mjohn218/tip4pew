CC =		mpic++

CCOPT =		-O3 -Wall #-showme 
#CCOPT =		-g

HOME =		.

DEFINE =	-DSGI_IRIX6_2

OBJS =		 Ensemble.o  


MYLIBS =	-lm

Tip4p_NVT: Tip4p_NVT.cpp $(OBJS) 
	$(CC) $(CCOPT) -o Tip4p_NVT Tip4p_NVT.cpp $(OBJS) $(MYLIBS)

nvt_nodiff: nvt_nodiff.cpp $(OBJS) 
	$(CC) $(CCOPT) -o nvt_nodiff nvt_nodiff.cpp $(OBJS) $(MYLIBS)

spce: spce.cpp $(OBJS) 
	$(CC) $(CCOPT) -o spce spce.cpp $(OBJS) $(MYLIBS)

clean:
	rm -f *.o Tip4p_NVT core 








