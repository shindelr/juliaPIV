"""
Top layer
"""

from .julia_compiler import compile_juliapiv
from .pivpipe.pivpipe import pivpipe_main, load_config
from .batcher.batcher import out_dir_setup, generate_txt_files
from .utils import batch_n_nproc_logic, build_dir_structure
import click
import os
import logging
from pathlib import Path
import yaml

logging.basicConfig(level=logging.INFO, format="%(asctime)s -- %(levelname)s -- %(message)s")

@click.group()
def juliaPIV():
    """
    Main entry for juliaPIV
    """
    pass

@juliaPIV.command()
def init():
    """
    Initialize the PIVPipelineUtility project. Compiling it into a full executable.
    """ 
    print("This may take a while...")
    compile_juliapiv()

@juliaPIV.command()
@click.option('--config', '-c', type=click.Path(exists=True), help="YAML config file to load defaults.")
def pipeline(config):
    """
    Run PIV pipeline.
    """
    settings = load_config(config)
    num_images = len(os.listdir(settings['input']))
    nproc, num_batches = batch_n_nproc_logic(settings['N'], settings["NPROC"], num_images)
    settings['NPROC'] = nproc

    logging.info("Building directory structure...")
    build_dir_structure(settings["output"])

    logging.info("Batching...")
    batches_dir = str(Path(settings["output"]).resolve() / "piv_batches")
    mat_dir = str(Path(settings["output"]).resolve() / "piv_mat_out")
    out_dir_setup(batches_dir)
    generate_txt_files(settings["input"], batches_dir, num_batches=num_batches)

    settings["input"], settings["output"] = batches_dir, mat_dir
    pivpipe_main(settings)
