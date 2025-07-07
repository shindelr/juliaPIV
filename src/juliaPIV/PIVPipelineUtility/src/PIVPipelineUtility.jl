module PIVPipelineUtility

using FileIO
using Images
using Statistics
using MAT
using ImageTransformations
include("./main.jl")

"""
    get_raw_images(path::String)::Vector{String}

    Read in a list of image names from a given text file. Note that if 
    the text file is not in the same directory as the images, the path to 
    the images must be prepended to the image names.   

    Arguments:
        - `path::String`: Relative path to the text file containing image names.

    Returns:
        - `Vector{String}`: Vector of image names.
"""
function get_raw_images(path::String, N::Int32)::Vector{String}
    files::Vector{String} = readlines(path)
    @assert length(files) % N == 0 "Number of files in batch ($(length(files))) should be divisible by N ($(N))" 
    return files
end

"""
    statistics_of_piv_groups(piv_results, N::Int32)

    Compute the statistics of PIV results for variable groups of
    images. After partitioning the results of the PIV algorithm into
    subgroups, each group's data is averaged and the standard deviation
    calculated.

    Arguments:
        - `piv_results`: An array containing each of the piv results 
            coming back from the main PIV algorithm.
    Returns:
        - `Tuple{Tuple{}, Tuple{}, Matrix{Int32}}`: A tuple of 
            tuples containing the averages, standard deviations, and a 
            sum of points that were NaN prior to averaging for each 
            pair of images.
"""
function statistics_of_piv_groups_memlite(piv_results, N::Int32)
    # Unpack results
    us, vs, = Vector{Matrix{Float32}}(), Vector{Matrix{Float32}}()
    xs, ys = Vector{Matrix{Float32}}(), Vector{Matrix{Float32}}()
    for result in piv_results
        push!(xs, result[1][1])
        push!(ys, result[1][2])
        push!(us, result[2][1])
        push!(vs, result[2][2])
    end
    
    nans = reduce(+, [isnan.(u) for u in us])
    u_avs = nan_mean(us)
    u_stds = nan_std(us)
    v_avs = nan_mean(vs)
    v_stds = nan_std(vs)
    npts = nans

    return ((xs[1], ys[1]), (u_avs, v_avs), (u_stds, v_stds), npts)
end

"""
    nan_mean(arr::Vector{Matrix{Float32}})::Matrix{Float32}

    Compute the mean of a vector of matrices, ignoring NaN values. 

    Arguments:
        - `arr::Vector{Matrix{Float32}}`: Vector of matrices containing 
            the data to be averaged. May contain NaN values.
    Returns:
        - `Matrix{Float32}`: A single matrix containing the mean of all the 
            input matrices disregarding NaN values.
"""
function nan_mean(arr::Vector{Matrix{Float32}})::Matrix{Float32}
    # Preallocate
    mean_matrix = Matrix{Float32}(undef, size(arr[1]))
    for i in 1:size(arr[1], 1)
        for j in 1:size(arr[1], 2)
            # If all values are NaN, mean_val = NaN
            mean_val = NaN
            temp_means = Vector{Float32}()
            for m in arr
                if !isnan(m[i, j])  # Only want the mean of non-nans
                    push!(temp_means, m[i, j])
                end
            end
            if length(temp_means) > 0
                mean_val = mean(temp_means)
            end
            mean_matrix[i, j] = mean_val
        end
    end
    @assert all(isnan.(mean_matrix)) == false "All values in mean_matrix should not be NaN"
    return mean_matrix
end

"""
    nan_std(arr::Vector{Matrix{Float32}})::Matrix{Float32}

    Compute the standard deviation of a vector of matrices, ignoring 
    NaN values. NOTE: If there is only one value in the vector, the
    standard deviation is undefined and set to 0.0. If all values in
    the matrix return NaN, the function will halt and throw an error.

    Arguments:
        - `arr::Vector{Matrix{Float32}}`: Vector of matrices containing 
            the data to be averaged. May contain NaN values.
    Returns:
        - `Matrix{Float32}`: A single matrix containing the std of all the 
            input matrices disregarding NaN values.
"""
function nan_std(arr::Vector{Matrix{Float32}})::Matrix{Float32}
    std_matrix = Matrix{Float32}(undef, size(arr[1]))
    for i in 1:size(arr[1], 1)
        for j in 1:size(arr[1], 2)
            std_val = NaN
            temp_stds = Vector{Float32}()
            for m in arr
                if !isnan(m[i, j])
                    push!(temp_stds, m[i, j])
                end
            end
            if length(temp_stds) > 1
                std_val = std(temp_stds)
            elseif length(temp_stds) == 1
                std_val = 0.0  # If only one value, std = 0 bc no variance for one value
            end
            std_matrix[i, j] = std_val
        end
    end
    @assert all(isnan.(std_matrix)) == false "All values in std_matrix should not be NaN"
    return std_matrix
end

"""
    Run PIV on a batch of images where N=2.
"""
function paired_piv(N::T, final_win_size::T, ol::Float32, out_dir::String, 
                    images::Vector{String}, crop_factor::Tuple{T,T,T,T}, 
                    downsample_factor::Float32, save_images::Bool) where {T}
    for i in 1:2:length(images) - 1
        name = replace(basename(images[i]), ".jpg" => "")
        img1 = Gray.(load(abspath(images[i])))
        img2 = Gray.(load(abspath(images[i+1])))

        # Crop/Downsample
        img1 = img1[crop_factor[3]:crop_factor[4], crop_factor[1]:crop_factor[2]]
        img2 = img2[crop_factor[3]:crop_factor[4], crop_factor[1]:crop_factor[2]]
        if downsample_factor < 1.0
            img1 = imresize(img1, ratio=downsample_factor)
            img2 = imresize(img2, ratio=downsample_factor)
        end

        if save_images
            save_path_dir = joinpath(out_dir, "cropped-downsampled-images")
            if !isdir(save_path_dir)
                mkdir(save_path_dir)
            end
            FileIO.save(joinpath(save_path_dir, basename(images[i])), img1)
            FileIO.save(joinpath(save_path_dir, basename(images[i+1])), img2)
        end

        # PIV!!
        try
            raw_piv_results = main((img1, img2), Int32(final_win_size), Float32(ol))
            println("Building .mat file --> $name")
            pass_sizes = [raw_piv_results[3] raw_piv_results[3]]  # Just a formatting thing to match OG Matlab
            for (i, result) in enumerate(raw_piv_results)
                x = raw_piv_results[1][1]
                y = raw_piv_results[1][2]
                u = raw_piv_results[2][1]
                v = raw_piv_results[2][2]
                npts = isnan.(u)
                mat_dict = Dict(
                    "x" => x,
                    "y" => y,
                    "pass_sizes" => pass_sizes,
                    "overlap" => ol,
                    "method" => "multin",
                    "fn" => name,
                    "u" => u,
                    "v" => v,
                    "npts" => npts,
                )
                MAT.matwrite("$out_dir/$name.mat", mat_dict)
            end
        catch e
            if isa(e, DimensionMismatch)
                @warn "Skipping because of dimension mismatch: $name"
                continue
            else
                throw(e)
            end
        end
    end

end


"""
    Run PIV when N > 2.
"""
function grouped_piv_memlite(N::T, final_win_size::T, ol::Float32, out_dir::String, 
                    images::Vector{String}, crop_factor::Tuple{T,T,T,T}, 
                    downsample_factor::Float32, save_images::Bool) where {T}

    for i in 1:N:length(images) - 1
        name = replace(basename(images[i]), ".jpg" => "")

        # Load/Crop/Downsample subset size N of images
        img_subset = Vector{Matrix{Gray{N0f8}}}()
        for j in 1:N
            img = Gray.(load(images[i + j - 1]))
            img = img[crop_factor[3]:crop_factor[4], crop_factor[1]:crop_factor[2]]
            if downsample_factor < 1.0
                img = imresize(img, ratio=downsample_factor)
            end
            push!(img_subset, img)

            # Save images if desired
            if save_images
                save_path_dir = joinpath(out_dir, "cropped-downsampled-images")
                if !isdir(save_path_dir)
                    mkdir(save_path_dir)
                end
                FileIO.save(joinpath(save_path_dir, basename(images[i+j-1])), img)
            end
        end
        try
            # Preallocate results from PIV: [(x, y), (u, v), pass_sizes]
            raw_piv_results = Vector{Tuple{
                                Tuple{Matrix{T}, Matrix{T}}, 
                                Tuple{Matrix{T}, Matrix{T}}, 
                                Vector{Int32}
                                } where {T}}()
            for i in (1:length(img_subset) - 1)
                push!(raw_piv_results, main((img_subset[i], img_subset[i+1]), Int32(final_win_size), Float32(ol)))
            end

            # PIV stats
            ((x, y),
            (u_av, v_av),
            (u_std, v_std), 
            npts) = statistics_of_piv_groups_memlite(raw_piv_results, N)

            println("Building .mat file --> $name")
            pass_sizes = [raw_piv_results[1][3] raw_piv_results[1][3]]  # Just a formatting thing to match OG Matlab
            npts = isnan.(u_av)
            mat_dict = Dict(
                "x" => x,
                "y" => y,
                "pass_sizes" => pass_sizes,
                "overlap" => ol,
                "method" => "multin",
                "fn" => name,
                "u" => u_av,
                "v" => v_av,
                "ustd" => u_std,
                "vstd" => v_std,
                "npts" => npts,
            )
            MAT.matwrite("$out_dir/$name.mat", mat_dict)
        catch e
            if isa(e, DimensionMismatch)
                @warn "Skipping because of dimension mismatch: $name"
                continue
            else
                throw(e)
            end
        end
    end
end

"""
    io_main(N::T, crop_factor::Tuple{T,T,T,T}, final_win_size::T,
    ol::Float32, out_dir::String, in_path::String) where {T}

    Run a single batch of images through the Julia PIV algorithm.

    Arguments:
        - `N::Int32`: Number of images in each subgroup to run PIV on.
            Corresponds to the LiDAR scan rate.
        - `crop_factor::Int32`: Tuple of 4 integers representing the 
            cropping factor (left, right, top, bottom).
        - `final_win_size::Int32`: Final window size for PIV.
        - `ol::Float32`: Overlap percentage for PIV.
        - `out_dir::String`: Directory to write .mat files to.
        - `in_path::String`: Path to the directory containing the images.

    Returns:
        None
"""
function io_main(N::T, crop_factor::Tuple{T,T,T,T}, final_win_size::T,
    ol::Float32, out_dir::String, in_path::String, downsample_factor::Float32, 
    save_images::Bool) where {T}

    # Image pre-processing
    images = get_raw_images(in_path, N)
    if length(images) % N != 0
        error("\n\nNumber of images in directory ($(length(images))) not divisible by N ($N).\n\n")
    end

    if N == 2
        paired_piv(N, final_win_size, ol, out_dir, images, crop_factor, downsample_factor, save_images)
        return 
    elseif N > 2
        # grouped_piv(N, final_win_size, ol, out_dir, images, crop_factor, downsample_factor)
        grouped_piv_memlite(N, final_win_size, ol, out_dir, images, crop_factor, downsample_factor, save_images)
        return
    else
        error("N should be greater than 1")
    end
end

"""
    Declares a function wrapper for io_main which allows C style input. This
    reduces overhead and should speed spin-up of the compiled Julia package.

    .mat files will be written
    to the argued `out_dir` directory path. These .mat files will contain
    a variety of information detailed here:
        x: [255x299 double]
        y: [255x299 double]
        pass_sizes: [3x2 double]
        overlap: 0.5
            method: 'multin'
            fn: {list of jpg files}
                u: [255x299 double]
                v: [255x299 double]
            npts: [255x299 double]  # number of data points that weren't NaN 
                                    # prior to time-average
            uStd: [255x299 double]  # standard deviation of the N results
            vStd: [255x299 double]  # ditto

    Returns:
        - `Cint`: 0 if successful, 1 if unsuccessful.
"""
function io_main_wrapper(
    N::Int32, 
    left::Int32, right::Int32, top::Int32, bottom::Int32,    # Crop factors
    final_win_size::Int32,
    ol::Float32, 
    out_dir::String, 
    in_path::String,
    quiet::Int,
    downsample_factor::Float32,
    save_images::Int,
    )::Cint

    # Run PIV pipeline
    try 
        crop_factors = (left, right, top, bottom)
        save_images_bool = Bool(save_images)
        
        # if quiet == 1     # Shut down stdout 
        #     og_stdout = stdout
        #     redirect_stdout(devnull)
        # end
        io_main(N, crop_factors, final_win_size, ol, out_dir, in_path, downsample_factor, save_images_bool)

        # if quiet == 1     # Bring it back
        #     redirect_stdout(stdout)
        # end

        return 0
    catch e
        @error "Error in io_main: $e"
        return 1
    end
end

end # end module