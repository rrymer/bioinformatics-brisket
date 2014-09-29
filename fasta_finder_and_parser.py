# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 14:25:28 2014

@author: pcrzz

These are the coding tasks 1 and 2
"""

import argparse
from twisted.python.filepath import FilePath
import os
from Bio import SeqIO as seqio

class Data(object):
    """
    A data object
    """
    def __init__(self,path,extension):
            """
            initializes a data object
            """
            self.basedir = path
            self.set_data_type(extension)
            self.list_of_filepaths = self.find_data_files()
            self.sequence_objects = self.get_sequence_objects()
    
    def set_data_type(self,extension):
        """
        Sets the data type of the Data object.  Also used as the extension for
        finding files.
        """
        self.data_type = extension.strip('.')
    
    def find_data_files(self):
        """
        returns a list of file paths
        """
        temp_list = []
        for path_name_tuple in os.walk(self.basedir):
            dir_name = path_name_tuple[0]
            for file_base_name in path_name_tuple[2]:
                if file_base_name.endswith('.' + self.data_type):
                    temp_list.append(FilePath(os.path.join(dir_name,file_base_name)))
        return temp_list
    
    def get_basedir(self):
        """
        """
        return self.basedir
    
    def get_sequence_objects(self):
        """
        Returns a dictionary with file_path objects as keys and sequence records
        as values.
        """
        sequence_records = {}
        for file_path in self.list_of_filepaths:
            sequence_records[file_path] = []
            record_iterator = seqio.parse(file_path.open(),self.data_type)
            for record in record_iterator:
                sequence_records[file_path].append(record)
        return sequence_records
    
    def get_sequences(self):
        """
        Returns a list of sequence strings
        """
        sequences = []
        for seq_record_list in self.sequence_objects.values():#should be sequence object lists
            for seq_record in seq_record_list:
                sequences.append(seq_record.seq)
        return sequences
        
    def get_percentage_read_length(self,length):
        """
        Determines the percentage of reads of length, length, in each data file
        """
        assert isinstance(length,int), 'must supply an integer read length'
        ratios = {}
        for file_path,sequences in self.sequence_objects.iteritems():
            ratios[file_path.basename()] = round(100*len(filter(lambda a: len(a.seq) >= length, sequences))\
            /float(len(sequences)),2)
        return self.print_statistic_by_file(ratios)    
    
    def print_statistic_by_file(self,data):
        """
        Takes a dictionary of statistics by file, and prints that statistic by 
        file in a readable format
        """
        for name,value in data.iteritems():
            print name,':', value
    
    def get_sequence_frequencies(self):
        """
        Returns the frequency of each sequence in the Data object
        """
        sequence_count_dict = {}
        count_sequence_dict = {}
        sequence_records = self.get_sequences()
        for sequence in sequence_records:
            try: 
                sequence_count_dict[sequence._data] += 1
            except KeyError:
                sequence_count_dict[sequence._data] = 1
        for key,value in sequence_count_dict.iteritems():
            count = value
            try:  
                count_sequence_dict[count].append(key)
            except:
                count_sequence_dict[count] = [key]
        i = 0
        while i < 10:
            max_count = max(count_sequence_dict.keys())
            print max_count,count_sequence_dict[max_count][0]
            count_sequence_dict.pop(max_count)
            i += 1
        
    def __str__(self):
        print 'Data object with', len(self.list_of_filepaths), 'data files.'

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Usage: python task1_2.py path/to/data/file <optional args>.\nTakes sequence files of user defined format, and calculates statistics on those files, such as the frequency of occurrence of each sequence.  If you have problems or questions, email Richard.Rymer@gmail.com.")
    parser.add_argument("src_path", metavar="path", type=str,
        help="Base directory of file paths containing sequences")
    parser.add_argument("--format",metavar="format", type=str, help="Format of the sequences files to be examined")
    parser.add_argument("--task", metavar="task", type=str, required=True,
        help="Input a task to complete, available options are: task1 and task2")

    args = parser.parse_args()
    file_format = args.format
    data = Data(args.src_path,file_format)
    task = args.task
    if task == 'task1':
        print 'Determining the percentage of reads of length 30 or greater for all files'
        data.get_percentage_read_length(30)
    elif task == 'task2':
        print 'Determining the 10 most common sequences in all files'
        data.get_sequence_frequencies()
    else:
        print 'That is not a valid task selection.'
    