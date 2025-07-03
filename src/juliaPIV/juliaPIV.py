"""
Top layer
"""

from .julia_compiler import compile_juliapiv
from .pivpipe.pivpipe import run_pipe
from .batcher.batcher import batcher_cli
from .utils import batch_n_nproc_logic
import click

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
    print(batch_n_nproc_logic(4, 10, 12))
    # run_pipe('test')

