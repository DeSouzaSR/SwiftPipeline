#!/usr/bin/env python
#-*- coding:utf-8 -*-
from configobj import ConfigObj
import os, sys
import shutil

def create_file(file):    
    simulation_name = file
    simulation_name = simulation_name.split('.')[0]

    if os.path.isdir(simulation_name):
        shutil.rmtree(simulation_name)
        os.mkdir(simulation_name)
    else:
        os.mkdir(simulation_name)
    
    with open(file, 'r') as f:
        config_file = yaml.load(f)

    suffix_simu = config_file['suffix_simu']
#    clones_numbers = config_file['clones_numbers']


    for i in range(1,len(suffix_simu) + 1):
        os.makedirs(str(simulation_name)+ "/" + str(simulation_name) + "{:03d}".format(i))
        
if __name__ == '__main__':
    create_file(sys.argv[1])
    
