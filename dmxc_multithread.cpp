#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <thread>
/*demultiplexing reads from a fastq file
based on i7 and i5 with corresponding
sample names in a BCsplit file.
assumes one line sequence in the fastq, six bp i7 and six bp i5 at end of header line format i7+i5
Written to work with GT-seq protocolas 
and data managment at EFGL
To execute, in command line type: dmxc ./path/to/BCsplit_file.csv ./path/to/library.fq
or type: dmxc -bc ./path/to/BCsplit_file.csv -lib ./path/to/library.fq -t num_threads_to_use
-t is optional, default is max number of threads available
Tom Delomas 2019*/

//function that performs the demultiplexing
void dmx(std::vector<std::string> output_files, std::vector<std::string> barcodes, std::string lib_path){

	//open library fastq
	std::ifstream library (lib_path);	
	
	//create barcode map and open output files streams
	std::unordered_map<std::string, std::ofstream*> barcode_map;	//create map of barcodes to pointers to output files
	
	for (int i=0; i<output_files.size(); i++){
		barcode_map[barcodes[i]] = new std::ofstream (output_files[i], std::ios::app);	//key is i7+i5 value is output file stream
	}
	
	//demultiplex
	std::string header;	//defining variables, one for each line in a read
	std::string read;
	std::string plus;
	std::string qual;
	std::string i7i5_str;	//string to store i7 and i5
	int length; //length of header string
	while (std::getline(library, header)){
		length = header.length();	//get length to get coordinates for slicing
		//get i7+i5
		i7i5_str=header.substr(length - 13, 13);
		//get other lines
		std::getline(library, read);
		std::getline(library, plus);
		std::getline(library, qual);
		//determine if barcodes match sample
		if (barcode_map.count(i7i5_str) == 1){
			*barcode_map[i7i5_str] << header << "\n" << read << "\n" << plus << "\n" << qual << "\n";
		}
	}
			
	//close library
	library.close();
	
	//close output files
	for(auto& pair : barcode_map) {
		delete pair.second;
		pair.second = 0;
	}
	
}

//reads barcode file, checks inputs, and spawns threads as appropriate
int main(int argc, char* argv[]){
	//get and check input options
	std::string bc_path = "NULL";
	std::string lib_path = "NULL";
	int max_thread = std::thread::hardware_concurrency();
	
	bool flags_used = false;	//check if flags were used for cmd line input or not
	std::string x;
	for(int i = 1; i < argc; i++){
		x = argv[i]; //convert to string
		if(x == "-bc"){
			flags_used = true;
			break;
		}
	}
	
	//determine command line input
	if (flags_used){
		for(int i = 1; i < argc; i++){
			x = argv[i]; //convert to string
			if(x == "-bc"){bc_path = argv[i+1];}
			if(x == "-lib"){lib_path = argv[i+1];}
			if(x == "-t"){
				std::string thread_str = argv[i+1];
				max_thread = stoi(thread_str);}
		}
	}
	else if (argc == 3){
		bc_path = argv[1];
		lib_path = argv[2];
	}
	else{
		std::cout<<"Error: incorrect number of arguments supplied on the command line. Exiting.\n";
		return 0;
	}
	
	if(bc_path.substr(bc_path.length() - 3, 3) != "csv"){
		std::cout<<"Error: BCsplit file is not a .csv file. Exiting.\n";
		return 0;
	}
	
	if(lib_path.substr(lib_path.length() - 2, 2) != "fq" && lib_path.substr(lib_path.length() - 5, 5) != "fastq"){
		std::cout<<"Error: library file is not a .fq or .fastq file. Exiting.\n";
		return 0;
	}
	if(max_thread < 1){
		std::cout<<"Error: number of threads is less than 1. Exiting.\n";
		return 0;
	}

	//open bcsplit file
	std::ifstream barcode_file (bc_path); 
	//declare vectors containing matching barcode and file names
	std::vector<std::string> all_output_files;
	std::vector<std::string> all_barcodes;
	
	if (barcode_file.is_open()){	
		//add barcodes and file names from BCsplit file
		std::string bc_line;
		std::vector<std::string> bc_cells;
		int pos;
		std::getline(barcode_file, bc_line);	//skip header row
		while (std::getline(barcode_file, bc_line)){
			if (bc_line.substr(0,1) != ","){	//check to make sure line represents a sample and is not a blank row
				bc_cells = {"", "", "", "", "", "", "", ""};
				pos=0;
				for(int i=0; i < bc_line.length(); i++){
					if(bc_line[i] == ','){pos++;}
					else if (bc_line[i] == '\n' || bc_line[i] == '\r'){}	//avoid including newline characters in variables
					else{bc_cells[pos] = bc_cells[pos] + bc_line[i];}
				}
				all_output_files.push_back(bc_cells[1] + bc_cells[0] + ".fastq");
				all_barcodes.push_back(bc_cells[5] + '+' + bc_cells[7]);
			}
		}
	}
	else{
		//file is not open, throw error and exit
		std:: cout << "could not open specified barcode file. exiting.\n";
		return 0;
	}
	
	//check that library can be opened
	//open library fastq
	std::ifstream library (lib_path);	
	
	if (!library.is_open()){
		//file is not open, throw error and exit
		std:: cout << "could not open specified library file. exiting.\n";
		return 0;
	} else {
		library.close();
	}
	int num_samples = all_output_files.size();
	
	std::cout << "Now demultiplexing " << num_samples << " sample(s) using up to " << max_thread << " thread(s).\nWhile you are waiting, ponder whether the following is true or false: \"This sentence is false.\"\n";

	//determine samples per thread
	int base_num = num_samples / max_thread; //integer division
	int remain = num_samples % max_thread;
	int count_done = 0; //count samples distributed to threads
	
	//spawn threads
	//example std::unordered_map<std::string, std::ofstream*> barcode_map;
	
	std::vector<std::thread> thread_vector;
	std::vector<std::string>::iterator output_iter = all_output_files.begin();	//iterators for output files and barcodes
	std::vector<std::string>::iterator barcode_iter = all_barcodes.begin();
	std::vector<std::string> temp_output_files;	//vectors to assign slices of output files and barodes and then pass to dmx()
	std::vector<std::string> temp_barcodes;
	
	//spawn threads with one extra sample
	for (int i=0; i < remain; i++){
		temp_output_files.assign(output_iter, output_iter + base_num + 1);	//slice of output files
		temp_barcodes.assign(barcode_iter, barcode_iter + base_num + 1);	//slice of barcodes
		output_iter = output_iter + base_num + 1;	//adjust iterators
		barcode_iter = barcode_iter + base_num + 1;
		//spawn thread
		thread_vector.push_back(std::thread(dmx, temp_output_files, temp_barcodes, lib_path));
		count_done = count_done + base_num + 1;
	}

	//spawn threads with base number
	while (count_done < num_samples){
		temp_output_files.assign(output_iter, output_iter + base_num);	//slice of output files
		temp_barcodes.assign(barcode_iter, barcode_iter + base_num);	//slice of barcodes
		output_iter = output_iter + base_num;	//adjust iterators
		barcode_iter = barcode_iter + base_num;
		//spawn thread
		thread_vector.push_back(std::thread(dmx, temp_output_files, temp_barcodes, lib_path));
		count_done = count_done + base_num;
	}
	
	//join threads
	for (auto& th : thread_vector){
		th.join();
	}
	
	std::cout<<"Demultiplexing completed.\n";

	return 0;
}
