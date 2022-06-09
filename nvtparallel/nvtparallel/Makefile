CC =		mpic++

CCOPT =		-O3 -Wall
#CCOPT =		-g

HOME =		.

DEFINE =	-DSGI_IRIX6_2

OBJS =		 Ensemble.o  


MYLIBS =	-lm

Tip4p_NVT: Tip4p_NVT.cpp $(OBJS) 
	$(CC) $(CCOPT) -o Tip4p_NVT Tip4p_NVT.cpp $(OBJS) $(MYLIBS)

nvt_diff2: nvt_diff2.cpp $(OBJS) 
	$(CC) $(CCOPT) -o nvt_diff2 nvt_diff2.cpp $(OBJS) $(MYLIBS)

nvt_diff: nvt_diff.cpp $(OBJS) 
	$(CC) $(CCOPT) -o nvt_diff nvt_diff.cpp $(OBJS) $(MYLIBS)

nvt_gr: nvt_gr.cpp $(OBJS) 
	$(CC) $(CCOPT) -o nvt_gr nvt_gr.cpp $(OBJS) $(MYLIBS)

pressurewrite: pressurewrite.cpp $(OBJS) 
	$(CC) $(CCOPT) -o pressurewrite pressurewrite.cpp $(OBJS) $(MYLIBS)

nvt_restart: nvt_restart.cpp $(OBJS) 
	$(CC) $(CCOPT) -o nvt_restart nvt_restart.cpp $(OBJS) $(MYLIBS)

Ensemble.o: Ensemble.cpp
	g++ -O3 -c Ensemble.cpp


clean:
	rm -f *.o Tip4p_NVT core 








