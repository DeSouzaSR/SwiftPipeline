#!/usr/bin/env python
#-*- coding:utf-8 -*-
import yaml
import os, sys 



def create_file():
    with open(sys.argv[1], 'r') as f:
        config_file = yaml.load(f)
    
    simulation_name = config_file['simulation_name']
    
    os.mkdir(simulation_name)
    
    
    
if __name__ == '__main__':
    create_file()
    
