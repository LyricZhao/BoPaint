APP = main

OBJECTS = main.o graph.o
OUTPUTS = output/*

CXXC = g++
CXXFLAGS = -std=c++11 -O2 -I/usr/local/opt/nlohmann_json/include
LINKFLAGS = -lopencv_core.3.4.3 -lopencv_highgui.3.4.3 -lopencv_imgproc.3.4.3 -lopencv_ml.3.4.3 -lopencv_imgcodecs.3.4.3

DEL = rm -rf

default:
	make clean -s
	make $(APP) -s

run:
	make $(APP) -s
	./$(APP) input/input.json output/output.png

$(APP): $(OBJECTS) Makefile
	echo 'Linking: $(APP)' && \
	$(CXXC) $(LINKFLAGS) $(OBJECTS) -o $(APP)

$(APP).o: $(APP).cpp Makefile
	echo 'Compiling main: $*.o' && \
	$(CXXC) $(CXXFLAGS) -c $*.cpp -o $*.o

%.o: %.cpp %.h Makefile
	echo 'Compiling components: $*.o' && \
	$(CXXC) $(CXXFLAGS) -c $*.cpp -o $*.o

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