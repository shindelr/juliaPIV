"""
Top layer
"""

from .julia_compiler import compile_juliapiv
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