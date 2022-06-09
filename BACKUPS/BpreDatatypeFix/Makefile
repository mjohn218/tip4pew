CC =		mpicxx

CCOPT =		-O3 -Wall
#CCOPT =		-g

HOME =		.

DEFINE =	-DSGI_IRIX6_2

OBJS =		 Ensemble.o  


MYLIBS =	-lm

Tip4p: Tip4p.cpp $(OBJS) 
	$(CC) $(CCOPT) -o Tip4p Tip4p.cpp $(OBJS) $(MYLIBS)

clean:
	rm -f *.o Tip4p core 








