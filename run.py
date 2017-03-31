#!/usr/bin/env python
#-*- coding:utf-8 -*-
import yaml
import system 

def create_file():
    with open('config.yaml', 'r') as f:
        config_file = yaml.load(f)
    
    simulation_name = config_file['simulation_name']
    
    system.mkdir(simulation_name)
    
    
    
if __name__ == '__main__':
    create_file()
    
