# Julia Particle Image Velocimetry

Welcome to the CEOAS Coastal Imaging Lab's RGB video PIV arm of the ROXSI processing pipeline! In this repository you can find a variety of directories related to PIV. This particular project was initiated by Prof. Wilson's desire to move over to an open-source programming platform for data processing. Julia checks all the boxes with a strong scientific programming community, JIT compilation leading to reasonably fast execution, voluntary type safety, and smooth read/write-ability.

By switching this process over to Julia, the execution time of processing a single pair of images improved by roughly 50%. Some improvements were made by tuning the overall algorithm and trimming unnecessary code. However, the majority of speed-up came from Julia's JIT compilation, typing, and some memory management since PIV is quite dependent on large arrays. 
