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
from tqdm import tqdm
from natsort import natsorted


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('indir')
@click.option('-o', '--out-dir', type=str, default='.',
              help='Directory to output .txt files')
@click.option('-n', '--n-files-per', type=int,
              help='Desired number of files per .txt output file, use this OR -b not both')
@click.option('-b', '--num-batches', type=int,
              help='Number of batches desired, use this OR -n not both')
def batcher_cli(indir, n_files_per, out_dir, num_batches):
    """
    A script to partition a directory of files into many directories with the OG
    files split between them. 
    """
    abs_in = os.path.abspath(indir)
    abs_out = os.path.abspath(out_dir)
    out_dir_setup(abs_out)
    generate_txt_files(abs_in, abs_out, n_files_per, num_batches)

def out_dir_setup(abs_out: str) -> None:
    """
    Ensure the out directory exists, creates it if it doesn't.
    """
    if not os.path.isdir(abs_out):
        print(f"Building a directory at {abs_out}")
        os.makedirs(abs_out)
        return
    if len(os.listdir(abs_out)) != 0:
        print("WARNING: piv_batches directory is not empty, exiting.")
        sys.exit()

def generate_txt_files(abs_in: str, abs_out: str, 
                       n_files_per: int = None,
                       num_batches: int = None
                       ) -> None:
    """
    Generate a number of .txt files based on how many individual files need to
    be parsed and how many files the user wants in each .txt file. Absolute paths
    to each file are built. 

    f-in-txt = f-total // n-files-per 
    """
    sorted_files = natsorted(os.listdir(abs_in))
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
