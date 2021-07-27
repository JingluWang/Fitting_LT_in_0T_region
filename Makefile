CXX = g++
ROOT_PATH = `root-config --cflags --libs ` -lMinuit
Fitting.exe : Fitting.C
	$(CXX) $(ROOT_PATH) $^ -o $@
