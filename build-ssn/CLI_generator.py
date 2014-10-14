# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 12:15:30 2014

@author: Richard Rymer

contains a class for the creation of simple and flexible command line calling
functions, hiding the necessary subprocess calls.
"""
import subprocess as sp

class CLI(object):
    def __init__(self,program,INPUT=None,OUTPUT=None,ERR=None):
        """
        Instantiates a CLI object.
        """
        self.input = INPUT
        self.output = OUTPUT
        self.err = ERR
        self.add_program(program)
    
    def add_program(self,program):
        """
        adds a program command line interface to the CLI object
        """
        setattr(self,program,self._generate_call_function(program))
    
    def add_program_options(self,program,options):
        """
        adds arbitary command line options to the program
        """
        getattr(self,program).options = options
    
    def _generate_call_function(self,program):
        """
        Returns a call function for the program
        """
        def call():
            print 'calling ' + program
            options = []
            for opt,val in getattr(self,program).options.items():
                options.append(opt)
                options.append(val)
            p = sp.Popen([program] + options,stdout=sp.PIPE)
            self.output = p.stdout
            return p
        return call

def example():
    cmd = CLI('makeblastdb')
    cmd.add_program_options('makeblastdb',{'-in':None,'-out':None,'-dbtype':None})
    cmd.makeblastdb.options['-in'] = '/tmp/bacteria.99.1.genomic.fna'
    cmd.makeblastdb.options['-out'] = '/tmp/test'
    cmd.makeblastdb.options['-dbtype'] = 'nucl'
    cmd.makeblastdb()
    for line in iter(cmd.output):
        print line
