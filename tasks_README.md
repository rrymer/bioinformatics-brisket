```
coding task scripts
task1_2.py
task3.py

VERSION
NA

SUMMARY
Two scripts completing the coding tasks

INSTALLATION
Dependencies: 
Python 2.7.1
	os
	sys
argparse 1.2.1
biopython 1.6.2
zope.interface 4.1.1 (twisted)


USAGES
Run from a terminal window with the following syntax:
task1_2.py:

python task1_2.py path/to/base-directory <additional args>

Takes sequence files of user defined format, and calculates statistics on those files, 
such as the frequency of occurrence of each sequence. If you have problems or questions, email
Richard.Rymer@gmail.com.

no output files

COMMAND LINE OPTIONS
  -h, --help            show this help message and exit
  --format format  Format of the sequences files to be examined
  --task task      Input a task to complete, available options are: task1 and
                   task2

task3.py:
Usage: python task3.py path/to/data/file --annotations=<path> --out_path=<path>. 

Takes a tab-delimited file of positions to annotate, and a gtf formatted annotation file,
and attemptes to annotate the input positions based on the data in the
annotation file. If you have problems or questions, email
Richard.Rymer@gmail.com.

outputs annotated_positions.txt containing the annotated positions to the out_path directory.

COMMAND LINE OPTIONS
  -h, --help            show this help message and exit
  --annotations ANNOTATIONS
                        define path to gtf formatted annotation file
  --out_path out_path   Define path for saving output files

DEVELOPMENT
NA

KNOWN ISSUES
NA

LEGAL
To do
```
