"""
Utility to compile the PIVPipelineUtility.jl and main.jl scripts.
"""

import subprocess
import os
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s -- %(levelname)s -- %(message)s")
logging.basicConfig(level=logging.ERROR, format="%(asctime)s -- %(levelname)s -- %(message)s")

def compile_juliapiv():
    """
    Compile juliaPIV on the target OS.
    """
    root_dir = Path(__file__).resolve().parent   # juliaPIV/src/juliaPIV/
    pivpipelineutility_path = os.path.join(root_dir, "PIVPipelineUtility")
    precompile_dir_path = os.path.join(pivpipelineutility_path, "src/precompile")
    precompile_script_path = os.path.join(precompile_dir_path, "precompile.jl")

    abs_paths_in_precompiletxt(precompile_dir_path)

    cmmd = [
        'julia',
        f'--project={pivpipelineutility_path}',
        '-e',
        f"""
        using Pkg;
        Pkg.add("PackageCompiler");
        using PackageCompiler;
        Pkg.instantiate();
        create_sysimage([:PIVPipelineUtility], sysimage_path="{os.path.join(pivpipelineutility_path, "pivbuildsysimage.so")}", precompile_execution_file="{precompile_script_path}", project="{pivpipelineutility_path}");
        """
    ]
    subprocess.run(cmmd)

    print(f"Build location: {os.path.join(pivpipelineutility_path, 'pivbuildsysimage.so')}")

def abs_paths_in_precompiletxt(precompile_dir):
    """
    Give each image absolute paths.
    """
    precompiletxt = Path(precompile_dir) / "in/batch/precompile.txt"
    jpg_dir = Path(precompile_dir) / "in/data"
    jpg_names = os.listdir(jpg_dir)
    jpgs_abs = [f"{str(jpg_dir.resolve() / jpg)}\n" for jpg in jpg_names]
    with open(precompiletxt, 'w') as f:
        f.writelines(jpgs_abs)
    