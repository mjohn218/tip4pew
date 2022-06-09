CC =		mpic++

CCOPT =		-O3 -Wall
#CCOPT =		-g

HOME =		.

DEFINE =	-DSGI_IRIX6_2

OBJS =		 md_timer.o Ensemble.o  


MYLIBS =	-lm

Tip4p: Tip4p.cpp $(OBJS) 
	$(CC) $(CCOPT) -o Tip4p Tip4p.cpp $(OBJS) $(MYLIBS)

Timetest: TimeTest.cpp $(OBJS) 
	$(CC) $(CCOPT) -o TimeTest TimeTest.cpp $(OBJS) $(MYLIBS)

Tip4pnewk: Tip4pnewk.cpp $(OBJS) 
	$(CC) $(CCOPT) -o Tip4pnewk Tip4pnewk.cpp $(OBJS) $(MYLIBS)

clean:
	rm -f *.o Tip4p core 








