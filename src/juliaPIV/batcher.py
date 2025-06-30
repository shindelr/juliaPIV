#!/Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11
"""
Author: Robin Shindelman
Date: 2025-04-09
A script to partition a directory of files into many directories with the OG
files split between them. This script was originally built to facilitate the
OSU ROXSI SVS pipeline.
"""

import os
import sys
import click
import re
from tqdm import tqdm

@click.command()
@click.argument('indir')
@click.option('-o', '--out-dir', type=str, default='.',
              help='Directory to output .txt files')
@click.option('-n', '--n-files-per', type=int,
              help='Desired number of files per .txt output file, use this OR -b not both')
@click.option('-b', '--num-batches', type=int,
              help='Number of batches desired, use this OR -n not both')
@click.option('-d', '--drone-sort', is_flag=True,
              help='Toggle on to accurately sort jpgs from Drone Footage')
def cli(indir, n_files_per, out_dir, num_batches, drone_sort):
    """
    A script to partition a directory of files into many directories with the OG
    files split between them. 
    """
    abs_in = os.path.abspath(indir)
    abs_out = os.path.abspath(out_dir)
    out_dir_setup(abs_out)

    if drone_sort:  # Check on file sort style
        generate_txt_files(abs_in, abs_out, n_files_per, num_batches, int_sort=False)
    else:
        generate_txt_files(abs_in, abs_out, n_files_per, num_batches, int_sort=True)


def out_dir_setup(abs_out: str) -> None:
    """
    Ensure the out directory exists, creates it if it doesn't.
    """
    if not os.path.isdir(abs_out):
        print(f"Building a directory at {abs_out}")
        os.makedirs(abs_out)
        return
    if len(os.listdir(abs_out)) != 0:
        print("WARNING: This directory is not empty, exiting.")
        sys.exit()



def generate_txt_files(abs_in: str, abs_out: str, 
                       n_files_per: int = None,
                       num_batches: int = None,
                       int_sort: str = True) -> None:
    """
    Generate a number of .txt files based on how many individual files need to
    be parsed and how many files the user wants in each .txt file. Absolute paths
    to each file are built. 

    f-in-txt = f-total // n-files-per 
    """
    key = int_sort_key if int_sort else drone_sort_key
    sorted_files = sorted(os.listdir(abs_in), key=key)
    files = [os.path.join(abs_in, f) for f in sorted_files]   # Absolute after sort
    num_f = len(files)

    if num_batches:
        f_per_batch = num_f // num_batches
        leftover = num_f % num_batches
        click.echo(f"Batches: {num_batches}; Files Per Batch: {f_per_batch}; Leftover: {leftover}; From {num_f} Files")

        with tqdm(total=num_f, desc="Batching") as pbar:
            chunk = 1
            i = 0
            while i < num_f and chunk <= num_batches:
                batch_fp = os.path.join(abs_out, f'batch{chunk}.txt')
                with open(batch_fp, 'a') as f:
                    f.write(f'{files[i]}\n')

                i += 1
                pbar.update(1)
                if i == (f_per_batch * chunk):  # Makes a new batch every n files
                    chunk += 1

    elif n_files_per:   # TODO: make this available
        click.echo("Dividing by desired number of files per batch")
        return
    else:
        click.echo("No batch number or file number defined. Aborting!")
    return


def drone_sort_key(f: str) -> int | float:
    """
    A sorting key for jpgs from ROXSI 2023. 
    Those files use the convention:
        A038_C002_1020WJ.015169.jpg
        A038_C002_1020WJ.015170.jpg
        A038_C002_1020WJ.015171.jpg
        A038_C002_1020WJ.015172.jpg
    
    REGEX: r'\.(\d+)(?=\.jpg$)'

    Matches on the ".numbers.jpg" group.
    """
    match = re.search(r'\.(\d+)(?=\.jpg$)', f)
    return int(match.group(1)) if match else float('inf')


def int_sort_key(f: str) -> int | float:
    """
    A sorting key for files beginning with a group of integers.
    Those fiels use the convention:
        011997_1738857742093682704_252_1190.jpg
        011998_1738857742106184504_252_1190.jpg
        011999_1738857742118686304_252_1190.jpg
        012000_1738857742131188096_252_1190.jpg
    
    REGEX: r'\d+'

    Simply matches on groups of numbers then takes the first match. If there are
    no digits in the filename, it is sorted to the bottom.
    """
    match = re.match(r'\d+', f)
    return int(match.group(0)) if match else float('inf')   # Matches on first group
