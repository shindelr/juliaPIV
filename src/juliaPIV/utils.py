"""
General utilities for the PIV pipeline.
"""

from multiprocessing import cpu_count

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