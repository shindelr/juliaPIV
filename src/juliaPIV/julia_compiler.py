"""
Utility to compile the PIVPipelineUtility.jl and main.jl scripts.
"""

import subprocess
import os
import pathlib
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
        create_sysimage([:PIVPipelineUtility], sysimage_path="{os.path.join(pivpipelineutility_path, "pivbuildsysimage.so")}", precompile_execution_file="{precompile_script_path}", project="{pivpipelineutility_path}");
        """
    ]
    subprocess.run(cmmd)

    print(f"Build location: {os.path.join(pivpipelineutility_path, 'pivbuildsysimage.so')}")
