# Julia Particle Image Velocimetry

Welcome to the CEOAS Coastal Imaging Lab's RGB video PIV arm of the ROXSI processing pipeline! In this repository you can find a variety of directories related to PIV. This particular project was initiated by Prof. Wilson's desire to move over to an open-source programming platform for data processing. Julia checks all the boxes with a strong scientific programming community, JIT compilation leading to reasonably fast execution, voluntary type safety, and smooth read/write-ability.

By switching this process over to Julia, the execution time of processing a single pair of images improved by roughly 50%. Some improvements were made by tuning the overall algorithm and trimming unnecessary code. However, the majority of speed-up came from Julia's JIT compilation, typing, and some memory management since PIV is quite dependent on large arrays. 


# Steps to first run on CEOAS servers
1. Ensure Python and Julia are loaded. A python version greater than or equal to 3.11 is required. Currently, the exact Julia version of 1.11.5 is required. 

    `ml eb-sw Julia`  
    `ml eb-sw Python`

2. In the directory of your choosing, create a Python virtual environment.

    `python3 -m venv .venv`

    For more information regarding Python virtual environments and package management, check out these links.  

    [Official Python Docs](https://docs.python.org/3/library/venv.html)  
    [Real Python](https://realpython.com/python-virtual-environments-a-primer/)  
This venv is where juliaPIV will be installed using `pip`.

**Step 1 must be done each time Monarch is fired up. Step 2 is a one time setup.**

3. Activate the virtual environment each time you'd like to use juliaPIV.

    `source .venv/bin/activate`

4. If this is the first time juliaPIV is being installed on your User, `pip` can
be used to get it:  
`pip install juliaPIV`  
Check that it was installed correctly with:  
`pip list`.

    To upgrade versions in the future, simply use:  
`pip install --upgrade juliaPIV`

5. **First time installations or new versions that required changes to Julia source code will require the following command.** The `init` command compiles the Julia code into a .so file, greatly decreasing system overhead during parallel processing and decreasing reliance on Just-In-Time compilation for common functions and data types. This process can take quite some time, and it is highly recommended that the command is **not** included in any automation scripts.

    `juliaPIV init`

6. Set up your desired configuration for PIV using a YAML file in the following format:

    ```YAML
    # config.yaml
    N: 2                                # Number of images in set to be processed
    crop_factor: "1, 6144, 1, 3240"     # Crop factors for image
    final_win_size: 16                  # Final window size for PIV
    ol: 0.5                             # Window overlap during PIV

    # Absolute path to directory containing images to be processed
    input: /home/monarch0/ROXSI/2025/pantera/flightData/20250214/Flight3_1303/splitjpg

    # Absolute path to desired output folder
    output: /home/server/pi/homes/chumar/piv_env

    quiet: false                        # Supresses print statements
    downsample_factor: 0.5              # Downsampling, float between 0 and 1
    NPROC: 40                           # Maximum number of parallel processes to run
    save_images: false                  # Save downsampled and cropped images
    ```

7. Run PIV!

    `juliaPIV pipeline -c path/to/your/yaml`
