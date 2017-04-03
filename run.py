#!/usr/bin/env python
#-*- coding:utf-8 -*-
import yaml
import os, sys

def create_file():
    with open(sys.argv[1], 'r') as f:
        config_file = yaml.load(f)
    
    simulation_name = sys.argv[1]
    simulation_name = simulation_name.split('.')[0]
    os.mkdir(simulation_name)
    
    os.mkdir()
if __name__ == '__main__':
    create_file()
    
