"""
General utilities for the PIV pipeline.
"""

from multiprocessing import cpu_count
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s -- %(levelname)s -- %(message)s")
logging.basicConfig(level=logging.WARNING, format="%(asctime)s -- %(levelname)s -- %(message)s")

def batch_n_nproc_logic(n, nproc, num_images):
    """
    Identify the best number of images to place in each batch depending on the values
    of N and NPROC.
    """
    settings_str = "NPROC: {}; N: {}; Number of images: {}"
    end_str = "NPROC: {}, N: {}, Batches: {}, Images per batch: {}, Total images to proces: {}"

    logging.info(f"User settings = \n\n\t{settings_str.format(nproc, n, num_images)}\n")

    if n < 2:
        logging.warning(f"N should be 2 or greater for PIV, aborting.")
        exit(1)

    max_cpus = cpu_count() - 1
    if nproc > max_cpus:
        logging.warning(f"NPROC set higher than the number of available cores: {max_cpus}")
        logging.warning(f"NPROC changed from {nproc} to {max_cpus}")
        nproc = max_cpus
    
    num_batches = find_num_batches(nproc, num_images, n)
    while num_images >= n: 
        num_batches = find_num_batches(nproc, num_images, n)
        if num_batches > 0:
            logging.info(f"Optimal settings: \n\n\t{end_str.format(num_batches, n, num_batches, num_images // num_batches, num_batches * (num_images // num_batches))}\n")
            return num_batches, num_batches,
        num_images -= 1

    logging.error(f"N = {n}, this may be more than the number of images!\n")
    raise ValueError("Could not find valid configuration given the settings.")

def find_num_batches(nproc, num_images, n):
    """
    TODO: Doc
    """
    for num_batches in range(nproc, 0, -1):
        img_p_batch = num_images // num_batches
        if img_p_batch >= n and img_p_batch % n == 0:
            return num_batches
    return 0

def build_dir_structure(parent):
    """
    Build the output directory structure.
    """
    parent = Path(parent).resolve()
    piv_mat_out = parent / "piv_mat_out"
    piv_batches = parent / "piv_batches"
    
    if not parent.is_dir():
        parent.mkdir()
    if not piv_mat_out.is_dir():
        piv_mat_out.mkdir()
    if not piv_batches.is_dir():
        piv_batches.mkdir()
