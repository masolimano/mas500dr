#!/usr/bin/env python
import os
import argparse
import numpy as np
import ccdproc as cdp



def main():
    pwd = os.getcwd()
    print(f'Current directory: {pwd}')
    print('This is script is under development.')
    print('Right now it only prints this message.')
    print('Soon it will serve to automatically reduce data from the MAS500 telescope.')
