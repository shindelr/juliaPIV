"""
Utility to compile the PIVPipelineUtility.jl and main.jl scripts.
"""

import subprocess
from platform import system
import os
import pathlib
from .utils import set_fallback_path_darwin
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s -- %(levelname)s -- %(message)s")
logging.basicConfig(level=logging.ERROR, format="%(asctime)s -- %(levelname)s -- %(message)s")

def compile_juliapiv():
    """
    Compile juliaPIV on the target OS.
    """
    root_dir = pathlib.Path(__file__).resolve().parent   # juliaPIV/src/juliaPIV/
    pivpipelineutility_path = os.path.join(root_dir, "PIVPipelineUtility")
    precompile_dir_path = os.path.join(pivpipelineutility_path, "src/precompile")
    precompile_script_path = os.path.join(precompile_dir_path, "precompile.jl")

    cmmd = [
        'julia',
        f'--project={pivpipelineutility_path}',
        '-e',
        f"""
        using Pkg;
        using PackageCompiler;
        Pkg.instantiate();
        create_library("{pivpipelineutility_path}", "src/juliaPIV/piv_build"; lib_name="pivbuild", precompile_execution_file="{precompile_script_path}", force=true);
        """
    ]
    subprocess.run(cmmd)

    print(f"Build location: {os.path.join(root_dir, 'piv_build')}")

def set_dyld_fallback():
    """
    DEPRECATED
    Set the DYLD fallback path. TODO: Add functionality for new OS's
    """
    # Set DYLD_FALLBACK_LIBRARY_PATH
    if system() == "Darwin":
        set_fallback_path_darwin()
    else:
        logging.error("OS not yet supported by juliaPIV")
        exit(1)