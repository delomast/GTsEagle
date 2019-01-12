# GTsEagle
A new pipeline for analyzing GT-seq data

Currently in development. Only the demultiplexing program (dmxc) is developed enough to be made public in its current form.

dmxc

dmxc quickly and efficiently dempultiplexes fastq files into separate fastq files based on i7 (index1) and i5 (index2) barcodes.

dmxc_multithread can be compiled with g++' using the following command:  
g++ dmxc_multithread.cpp -o dmxc -pthread

The resulting dmxc requires two inputs with an optional third input that is the maximum number of threads to use. The two required inputs are the library and the BCsplit file. The library is the fastq file that contains reads from all the samples you want to demultiplex. The BCsplit file is a comma separated value file that contains several fields of information about each sample.

The BCsplit file must have eight fields and follow the below format with each sample being a new line (the header is REQUIRED):

SampleName,SampleType,SampleStatus,PlateID,i7_name,i7_sequence,i5_name,i5_sequence  
sample_1001,initial,normal,Tray30,i101,TCTCGA,3,CGGAAT  
sample_1002,rerun,normal,Tray30,i101,TCTCGA,4,GCCTCG  

All fields must be present, but only SampleName, SampleType, i7_sequence, and i5_sequence are used, so the other fields can simply be filled with placeholder values. dmxc will create a fastq file named 'SampleTypeSampleName.fastq' for each  line of the BCsplit file with reads matching the corresponding i7 and i5 sequences. Work is planned to allow the program to accept an alternative format of the BCsplit file that is simpler - the current format accepted is a file that is used in my current laboratory, and so it was the first format I developed the program to accept. The i7 and i5 sequences are the sequences expected IN THE LIBRARY FILE, so depending on the machine you use, the i5 will either be the original sequence or the reverse complement. Currently, dmxc only excepts barcodes that are 6bp long, and only functions with libraries whose barcodes are 6bp long each and are at the end of the fastq header line with the format AAAAAA+CCCCCC


To run dmxc, use one of the following syntax options:

dmxc /path/to/BCsplit_file.csv /path/to/library.fq  
  Note: using this syntax implies using as many threads as your system has cores
  
or

dmxc -bc /path/to/BCsplit_file.csv -lib /path/to/library.fq -t max_number_threads_to_use  
  Note: the -t is optional, if it is not specified, then this implies using as many threads as your system has cores
  
