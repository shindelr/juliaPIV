"""
General utilities for the PIV pipeline.
"""

from multiprocessing import cpu_count
import logging
import subprocess
import os
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s -- %(levelname)s -- %(message)s")

def batch_n_nproc_logic(n, nproc, num_images):
    """
    Identify the best number of images to place in each batch depending on the values
    of N and NPROC.
    """
    header_string = "NPROC: {}; N: {}; Number of batches: {}"

    num_cpus = cpu_count()
    if (num_images % n == 0) and (nproc <= (num_cpus - 1)):
        return nproc
    if (num_images % n == 0) and (nproc > (num_cpus -1)):
        print(f"Reduced NPROC to Maximum Available CPUs: {num_cpus - 1}")
        print(header_string.format(num_cpus-1, n, num_cpus-1))
        return num_cpus - 1

    # if (num_images % n != 0) and 

def set_fallback_path_darwin():
    """
    DEPRECATED
    Set the DYLD fallback library path for julia on a MacOS.
    """
    # Find libjulia-internal first
    find_expanded = os.path.expanduser("~/.julia/juliaup")
    find_cmmd = ["find", find_expanded, "-name", "libjulia-internal.dylib"]
    find_out = subprocess.run(find_cmmd, capture_output=True)
    if find_out.returncode != 0:
        logging.error("juliaup not found.")
        exit(1)
    fallback_lib_path = Path(find_out.stdout.decode('utf-8')).parent
    print(fallback_lib_path)

    # Export the julia path into ENV
    env_query = os.environ.get("DYLD_FALLBACK_LIBRARY_PATH", "").split(":")
    for path in env_query:
        if str(fallback_lib_path) not in path:  # Make sure it doesn't exist
            print(f'{"="*160}\n\n')
            print("\tYou will only see this message once per setup")
            print("\tRun this terminal command:\n")
            print(f"\t\texport DYLD_FALLBACK_LIBRARY_PATH={fallback_lib_path}:$DYLD_FALLBACK_LIBRARY_PATH\n\n")
            print(f'{"="*160}\n\n')
        return
