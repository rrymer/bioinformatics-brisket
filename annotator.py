# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 14:25:28 2014

@author: pcrzz

This is coding task 3
"""
import argparse
import sys
import pdb

class Data(object):
    def __init__(self,path):
            """
            initializes a data object
            """
            self.positions = []
            self.chromosome_ranges = {}
            self.annotated_positions = []
            with open(path,'r') as fin:
                for line in fin:
                    self.add_position(line.strip().split('\t'))
    
    def add_position(self,coord):
        """
        instantiates a Position object with coordinates (chromosome, nt position),
        and adds that object to a list of Position objects
        """
        self.add_position_to_chromosome_range(coord)
        self.positions.append(Position(coord))

    def add_position_to_chromosome_range(self,coord):
        """
        Takes a coordinate, and adds the position to a dict covering all positions
        for each chromosome in the data set
        """
        try:
            self.chromosome_ranges[coord[0]].append(int(coord[1]))
        except KeyError:
            self.chromosome_ranges[coord[0]] = [int((coord[1]))]
    
    def get_chromosome_range(self,chromosome):
        """
        Returns the min and max positions as a tuple for a chromosome in the data set
        """
        return min(self.chromosome_ranges[chromosome]),max(self.chromosome_ranges[chromosome])
    
    def get_all_positions(self):
        """
        Returns all of the Position objects associated with the Data object
        """
        return self.positions

class Position(object):
    def __init__(self,coords):
        """
        instantiates an unannotated Position object with attribute coords, 
        where coords are a list or tuple of form: (chromosome, position)
        """
        self.chromosome = coords[0]
        self.position = int(coords[1])
        self.annotated = False
    
    def set_annotated_flag(self,value):
        """
        Sets the annotated status of the Position object
        """
        self.annotated = value
    
    def add_annotation(self,annotation):
        """
        Takes a gene name, id, etc., and adds it as an annotation to the Position 
        object as an attribute, annotation
        """
        if annotation:
            self.annotation = annotation
            self.set_annotated_flag(True)
        else:
            self.set_annotated_flag(False)
    
    def get_chromosome(self):
        """
        Returns the chromosome associated with the position
        """
        return self.chromosome
    
    def get_position(self):
        """
        Returns the nucleotide position of the Position object
        """
        return self.position
    
    def get_annotation(self):
        """
        Returns the annotation of the Position
        """
        if self.annotated:
            return self.annotation
        else:
            return ' Not Annotated'
    
    def is_annotated(self):
        """
        Returns the annotation status of the Position
        """
        return self.annotated
        
class Chromosome(object):
    def __init__(self,name,gene_name,gene_range):
        """
        Instantiates a Chromosome object with the nucleotide position range data
        for one gene
        """
        self.name = name
        self.ranges_of_genes = {}
        self.add_gene(gene_name,gene_range)
        self.sorted = False
    
    def add_gene(self,name,gene_range):
        """
        Adds a gene and nucleotide range (tuple) to the Chromosome object
        """
        self.ranges_of_genes[(int(gene_range[0]),int(gene_range[1]))] = name
    
    def create_sorted_ranges(self):
        """
        Sorts the range keys, and stores them in a separate list object
        """
        self.sorted_ranges = sorted(self.get_all_ranges(), key=lambda tup: tup[1])
        self.sorted = True
    
    def get_all_ranges(self):
        """
        Returns all annotated gene ranges for the chromosome
        """
        if self.was_sorted():
            return self.sorted_ranges
        else:
            return self.ranges_of_genes.keys()
    
    def get_gene_of_range(self,range_tuple):
        """
        Given a gene range tuple, returns the name of the gene for that region
        """
        return self.ranges_of_genes[range_tuple]
    
    def was_sorted(self):
        """
        Returns whether or not the chromosome range list has been sorted
        """
        return self.sorted

class Annotations(object):
    def __init__(self,path):
        """
        Initializes an annotation file (gtf formatted) object
        """
        self.chromosomes = {}
        with open(path,'r') as fin:
            for line in fin:
                split_line = line.strip().split('\t')
                assert len(split_line) == 9
                try:
                    self.chromosomes[split_line[0]].add_gene(clean_up_annotation(split_line[-1]),
                                                                (split_line[3],split_line[4]))
                except KeyError:
                    self.chromosomes[split_line[0]] = Chromosome(split_line[0],clean_up_annotation(split_line[-1]),
                                                                (split_line[3],split_line[4]))
        self.sort_all_chromosomes()
        
    def sort_all_chromosomes(self):
        """
        Iterates over all of the chromosome objects in the annotations data set, and sorts the range lists
        """
        for chromosome in self.chromosomes.values():
            chromosome.create_sorted_ranges()
    
    def get_chromosome_annotations(self,chromosome):
        """
        takes a chromosome id as a string, and returns the genes annotated to that
        chromosome
        """
        if chromosome == 'chrM':
            return self.chromosomes['chrY']
        else:
            return self.chromosomes[chromosome]

def update_progress(num,max_progress,ratio):
    """
    generates a progress bar, based on: http://stackoverflow.com/questions/3160699/python-progress-bar
    """
    if num == 1 or 100*(float(num)/(max_progress))%1 <= ratio:
        pass
    else:
        return
    progress = float(num)/max_progress
    barLength = 100 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), round(progress*100,0), status)
    sys.stdout.write(text)
    sys.stdout.flush()

def clean_up_annotation(annotation):
    """
    Takes a raw annotation and returns just the gene name
    """
    return annotation.split(';')[-2].split(' ')[-1].strip('\"')

def binary_search(position,ranges):
    """
    Given a position and list of ranges, finds the relative locations of that
    position in the ranges list using a binary search algorithm
    """
    if len(ranges) <= 1:
        if position >= ranges[0][0] and position <= ranges[0][1]:
            return ranges[0]
        else:
            return None  
    elif position >= ranges[len(ranges)/2] and position <= [len(ranges)/2][1]:
        return ranges[len(ranges)/2]
    elif position <= ranges[len(ranges)/2][0]:
        ranges = ranges[:len(ranges)/2]
        return binary_search(position,ranges)
    elif position >= ranges[len(ranges)/2][1]:
        ranges = ranges[len(ranges)/2:]
        return binary_search(position,ranges)
    elif len(ranges) <= 3:
        for gene_range in ranges:
            if position >= gene_range[0] and position <= gene_range[1]:
                return gene_range
            else:
                continue
        return None
    

def find_gene_from_position(position,annotations,data):
    """
    Takes a Position object and an Annotations object, and searches through
    the annotations for a possible gene match, which is added to the Position
    as an annotation.
    """
    chromosome = annotations.get_chromosome_annotations(position.get_chromosome())
    assert chromosome.get_all_ranges()[0][0] < chromosome.get_all_ranges()[-1][1]
#    pdb.set_trace()
    if position.get_position() <= chromosome.get_all_ranges()[0][0] \
        or position.get_position() >= chromosome.get_all_ranges()[-1][1]:
           return None
    annotated_range = binary_search(position.get_position(),chromosome.get_all_ranges())
    if annotated_range:
        return chromosome.get_gene_of_range(annotated_range)
    else:
        possibles = chromosome.get_all_ranges()
        for annotated_range in possibles:
            if int(position.get_position()) >= int(annotated_range[0]) and int(position.get_position()) <= int(annotated_range[1]):
                return chromosome.get_gene_of_range(annotated_range)
            else:
                continue
        return None
            
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="python task3.py Usage: path/to/data/file --annotations=<path> --out_path=<path>.\nTakes a tab-delimited file of positions to annotate, and a gtf formatted annotation file, and attemptes to annotate the input positions based on the data in the annotation file. If you have problems or questions, email Richard.Rymer@gmail.com.")
    parser.add_argument("--annotations" ,dest="annotations", required=True,
                            help="define path to gtf formatted annotation file")
    parser.add_argument("src_path", metavar="path", type=str,
        help="Path to file containing the positions to annotate, only one file maybe processed at a time")
    parser.add_argument("--out_path", metavar="out_path", type=str,
        help="Define path for saving output files")

    args = parser.parse_args()
    data = Data(args.src_path)
    out_path = args.out_path
    annotations = Annotations(args.annotations)
    sys.setrecursionlimit(100000)
    i = 0
    max_positions = len(data.get_all_positions())
    for position in data.get_all_positions():
        position.add_annotation(find_gene_from_position(position,annotations,data))
        update_progress(i,max_positions-1,0.1)
        i += 1
    annotated = 0
    missed = 0
    with open(out_path + '/coordinates_annotated.txt','w') as fout:
        for position in data.get_all_positions():
            if position:
                if position.is_annotated():
                    annotated += 1
                else:
                    missed += 1
                fout.write('{}\t{}\t{}'.format(position.get_chromosome(),position.get_position(),position.get_annotation()) + '\n')
    print '{} total positions assessed, of which {} positions were annotated, {} positions were missed'.format(max_positions,annotated,missed)
