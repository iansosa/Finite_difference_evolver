MAIN?=FDE
INCLUDE_ENSEMBLE?=Ensemble_param

all: compile clean_o run

compile: ${MAIN} ${INCLUDE_ENSEMBLE}
	g++ ${MAIN}.o ${INCLUDE_ENSEMBLE}.o -o ${MAIN}_run -w -fopenmp 

${INCLUDE_ENSEMBLE}: $(INCLUDE_ENSEMBLE).cpp
	g++ -c $< -o $@.o -w -fopenmp

${MAIN}: $(MAIN).cpp 
	g++ -c $< -o $@.o -w -fopenmp

run:
	./$(MAIN)_run

clean_o:
	rm *.o

clean:
	rm ${MAIN}_run