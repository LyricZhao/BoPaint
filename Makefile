APP = main

OBJECTS = main.o
OUTPUTS = output.png

CXXC = g++
CXXFLAGS = -std=c++11 -O2
LINKFLAGS = -lopencv_core.3.4.3 -lopencv_highgui.3.4.3 -lopencv_imgproc.3.4.3 -lopencv_ml.3.4.3 -lopencv_imgcodecs.3.4.3

DEL = rm -rf

default:
	make clean -s
	make $(APP) -s

$(APP): $(OBJECTS) Makefile
	echo 'Linking: $(APP)' && \
	$(CXXC) $(LINKFLAGS) $(OBJECTS) -o $(APP)

%.o: %.cpp Makefile
	echo 'Compiling: $*.o' && \
	$(CXXC) $(CXXFLAGS) -c $*.cpp -o $*.o

%.o: %.cpp %.h Makefile
	echo 'Compiling: $*.o' && \
	$(CXXC) $(CXXFLAGS) -c $*.cpp -o $*.o

run:
	./$(APP)

clean:
	echo 'Cleaning all files ...' 
	make -s clean_object
	make -s clean_output

clean_object:
	echo 'Cleaning objects ...'
	-$(DEL) *.o
	-$(DEL) $(APP_NAME)

clean_output:
	echo 'Cleaning outputs ...'
	-$(DEL) $(OUTPUTS)	