#!/usr/bin/env python
#-*- coding:utf-8 -*-
import yaml
import os, sys

def create_file(file):    
    simulation_name = file
    simulation_name = simulation_name.split('.')[0]
    os.mkdir(simulation_name)
    
    with open(file, 'r') as f:
        config_file = yaml.load(f)

    suffix_simu = config_file['suffix_simu']
    clones_numbers = config_file['clones_numbers']
    
    for i in range(len(suffix_simu)):
        os.mkdirs(str(simulation_name) + str(simulation_name) + str(i))
        
if __name__ == '__main__':
    create_file(sys.argv[1])
    
