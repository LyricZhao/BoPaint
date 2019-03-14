/* 
 * Author: Chenggang Zhao, CST 75, Tsinghua University
 */

# define APP_VERSION "1.0"

# include "graph.h"

# include <iostream>
# include <string>

void printHelp() {
	std:: cout << "BOPAINT v" + std:: string(APP_VERSION) + " ;" << std:: endl;
	std:: cout << "Usage: bopaint [input json] [output filename] ;" << std:: endl;
	std:: cout << "Help for json format is in readme ;" << std:: endl;
	return;
}

int main(int argc, char** argv) {
	if (argc != 3) {
		printHelp();
	} else {
		std:: string input_file = std:: string(argv[1]);
		std:: string output_file = std:: string(argv[2]);
		Graph graph(input_file);
		graph.draw();
		graph.saveToFile(output_file);
	}
	return 0;
}