#include <vector>
#include <string>
#include <iostream>
#include <fstream>

// Function used to calculate the short range wake.
// a is the array that holds the time info of all particles,
// b is the array that holds 
void dump_to_file(std::string& path, std::vector<double>& data){
    std::filebuf fb;
	fb.open(path, std::ios::out);
	std::ostream os(&fb);
	for (int i = 0; i<data.size(); ++i) {
		os << data[i] << '\n';
	}
	fb.close();
};
