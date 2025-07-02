# Third party modules
using Statistics
using FFTW            # Fast Fourier Transforms library built on C
using Images          # Basic image processing library
# using FileIO          # I/O library
# using DelimitedFiles  # Write matrices to CSV
using Interpolations
# using Plots
using Luxor            # For creating inpolygon() functionality
# include("./threaded_1stpass.jl")

# MAIN
"""
    main(image_pair::Tuple{Matrix{T}, Matrix{T}}) where {T}

    Main Entry
    Minimal PIV calculation for testing and development.
    Run two images through the PIV process, eventually ending in the expected PIV
    plots.
    This should be the only public facing function in this library.
    
    Parameters
    ----------

    Returns
    ----------

"""
function main(image_pair::Tuple{Matrix{T},Matrix{T}}, 
                final_win_size::Int32, 
                ol::Float32) where {T}
    # Convert the images to matrices of floats
    A = convert(Matrix{Float32}, image_pair[1])
    B = convert(Matrix{Float32}, image_pair[2])

    pass_sizes = zeros(Int32, 4)
    try 
        log2pivwin::Int32 = Int32(log2(final_win_size))
        pass_sizes::Vector{Int32} = 2 .^ (Int32(6):-1:log2pivwin)
        push!(pass_sizes, final_win_size) # Duplicate final element
    catch 
        InexactError
        error("final_win_size must be factor of 2")
    end

    # other input params for piv
    dt::Int32 = 1
    overlap::Float32 = ol
    validvec::Int32 = 3
    x, y, u, v, SnR, Pkh = multipassx(A, B, pass_sizes, dt, overlap, validvec)

    # Reject data with too-low signal-to-noise level
    snrthresh::Float32 = 1.3
    u .= ifelse.(SnR .< snrthresh, NaN, u)
    v .= ifelse.(SnR .< snrthresh, NaN, v)

    # Reject data with too-low correlation peak height
    pkhthresh::Float32 = 0.3
    u .= ifelse.(Pkh .< pkhthresh, NaN, u)
    v .= ifelse.(Pkh .< pkhthresh, NaN, v)

    # Reject data that disagree strongly with their neighbors in a local window
    u, v = globfilt(u, v)

    return ((x, y), (u, v), pass_sizes)

    # # Plotting stuff
    # u_map = heatmap(u, 
    #                 title = "u [pixels/frame]", 
    #                 aspect_ratio = :equal, 
    #                 limits=(0, 200), 
    #                 xlimits=(0, 385))

    # v_map = heatmap(v, 
    #                 title = "v [pixels/frame]", 
    #                 aspect_ratio = :equal, 
    #                 ylimits=(0, 200), 
    #                 xlimits=(0, 385))
    # dbl_plot = plot(u_map, v_map, layout = (2, 1))
    # dest_dir = "/home/server/pi/homes/shindelr/2025-06-27-fligh1-error-testing"
    # png(dbl_plot, dest_dir * "037047-037048-out.png")
end

# PASS FUNCTIONS 
"""
    multipassx(A::Matrix{Float32}, B::Matrix{Float32}, wins::Vector{Int32}, Dt::Int32, 
                    overlap::Float32, sensit::Int32)

    Parameters
    ----------
    A: Matrix containing image data of first frame.
    B: Matrix containing image data of second frame.
    wins: 2D matrix of ints containing sub-window pixel sizes for each pass.
    Dt: Frame time step in seconds (int). Pass 1 for 'pixels per frame velocity'.
    overlap: Fraction of window overlap. Int.
    sensit: Threshold for vector validation. Int. 

    Returns
    ----------
    u: U-component of pixel velocity
    v: V-component of pixel velocity
    SnR: Ratio representing signal-to-noise.
    Pkh: Peak height for use in validation of vector field?
"""
function multipassx(A::Matrix{T}, B::Matrix{T}, wins::Vector{Int32}, Dt::Int32,
    overlap::Float32, sensit::Int32) where {T}
    sy, sx = Int32.(size(A))
    total_passes = Int32(length(wins))

    # Initial passes are for removing large-scale displacements. Initialize
    # displacements (datax,datay) to zero
    data_dim_1::Int32 = round(Int32, sy / (wins[1] * (1 - overlap)))
    data_dim_2::Int32 = round(Int32, sx / (wins[1] * (1 - overlap)))
    datax = zeros(Float32, (data_dim_1, data_dim_2))
    datay = copy(datax)

    for i in 1:total_passes-1
        # i = 1
        # println("Pass ", i, " of ", total_passes)
        x, y, datax, datay = firstpass(A, B, wins[i], overlap, datax, datay)
        # @time x, y, datax, datay = threaded_firstpass(A, B, wins[i], overlap, datax, datay)
        # display(datax)
        
        datax, datay = localfilt(x, y, datax, datay, sensit)
        datax, datay = linear_naninterp(datax, datay)

        datax = floor.(datax)
        datay = floor.(datay)

        if i != total_passes - 1
            Y, X, YI, XI = build_grids(overlap, sy, sx, wins[i], wins[end])
            # Y, X, YI, XI = build_grids_2(datax)
            datax = regular_interp(datax, X, Y, XI, YI)
            datay = regular_interp(datay, X, Y, XI, YI)

            # datax = make_nan_border(datax)
            # datay = make_nan_border(datay)

            datax, datay = linear_naninterp(datax, datay)
            datax = round.(datax)
            datay = round.(datay)
        end
    end

    # println("Final Pass")
    x, y, u, v, SnR, Pkh = finalpass(A, B, wins[end], overlap, datax, datay, Dt)
    return x, y, u, v, SnR, Pkh
    # return 0, 0, 0, 0, 0, 0
end

"""
### firstpass
    *Note: First pass is a misnomer for this function, as it's called N-1 times
    depending on the size of the given windows. Consider renaming*\n\n
    Set up matrix indices along the image frame according to the desired overlap and \
    window sizes. Calls **xcorrf2** which finally uses FFTs to calculate cross \
    correlation.
    \n**:params:**\n
    A: Matrix containing image data of first frame.\n
    B: Matrix containing image data of second frame.\n
    N: Vector representing sub-window pixel sizes for the pass.\n
    overlap: Fraction of window overlap. Int.\n
    idx: Matrix of same type as A, containing data displacement information.\n
    idy: Matrix of same type as A, containing data displacement information.\n
"""
function firstpass(A::Matrix{T}, B::Matrix{T}, N::Int32, overlap::Float32,
    idx::Matrix{T}, idy::Matrix{Float32}) where {T}
    M = N

    # Set up for FFT plans
    pad_matrix_a = pad_for_xcorr(A[1:M, 1:N])
    pad_matrix_b = pad_for_xcorr(B[1:M, 1:N])
    P = plan_fft(pad_matrix_a; flags=FFTW.ESTIMATE)
    Pi = plan_ifft(pad_matrix_a; flags=FFTW.ESTIMATE)

    # Initializing matrices
    sy, sx = size(A)
    xx_dim1 = ceil(Int32, ((size(A, 1) - N) / ((1 - overlap) * N))) + 1
    xx_dim2 = ceil(Int32, ((size(A, 2) - M) / ((1 - overlap) * M))) + 1
    datax = zeros(eltype(A), (xx_dim1, xx_dim2))
    datay = zeros(eltype(A), (xx_dim1, xx_dim2))
    xx = zeros(eltype(A), (xx_dim1, xx_dim2))
    yy = zeros(eltype(A), (xx_dim1, xx_dim2))

    cj = 1
    for jj in 1:((1-overlap)*N):(sy-N+1)
        ci = 1
        for ii in 1:((1-overlap)*M):(sx-M+1)
            if isnan(idx[cj, ci])
                idx[cj, ci] = 0
            end

            if isnan(idy[cj, ci])
                idy[cj, ci] = 0
            end

            if (jj + idy[cj, ci]) < 1
                idy[cj, ci] = 1 - jj
            elseif (jj + idy[cj, ci]) > (sy - N + 1)
                idy[cj, ci] = sy - N + 1 - jj
            end

            if (ii + idx[cj, ci]) < 1
                idx[cj, ci] = 1 - ii
            elseif (ii + idx[cj, ci]) > (sx - M + 1)
                idx[cj, ci] = sx - M + 1 - ii
            end

            # Get windows
            C = A[floor(Int32, jj):floor(Int32, jj + N - 1),
                floor(Int32, ii):floor(Int32, ii + M - 1)]
            D = B[floor(Int32, jj + idy[cj, ci]):floor(Int32, jj + N - 1 + idy[cj, ci]),
                floor(Int32, ii + idx[cj, ci]):floor(Int32, ii + M - 1 + idx[cj, ci])]

            C = C .- mean(C)
            D = D .- mean(D)

            # Call xcorrf2, passing in the FFT plans
            R::Matrix{Float32} = xcorrf2(C, D, P, Pi, pad_matrix_a, pad_matrix_b)

            # Find position of maximal value of R
            max_coords = Vector{NTuple{2, Float32}}()
            subset = R[Int32(0.5 * N + 2):Int32(1.5 * N - 3), Int32(0.5 * M + 2):Int32(1.5 * M - 3)]
            fast_max!(max_coords, subset)

            # Adjust for subset positions
            max_coords = [(i[1] + Int32(0.5 * N + 1), 
                           i[2] + Int32(0.5 * M + 1))
                           for i in max_coords]

            # Handle a vector that has multiple maximum coordinates. Take the 
            # weighted average of the coordinates.
            if length(max_coords) > 1
                max_x1 = round(Int32, sum([c[2] * i for (i, c) in enumerate(max_coords)]) / sum([c[2] for c in max_coords]))
                max_y1 = round(Int32, sum([c[1] * i for (i, c) in enumerate(max_coords)]) / sum([c[1] for c in max_coords]))

            elseif isempty(max_coords)
                idx[cj, ci] = NaN
                idy[cj, ci] = NaN
                max_x1 = NaN
                max_y1 = NaN

            # Otherwise, unpack into max coordinates
            else
                max_y1, max_x1 = max_coords[1][1], max_coords[1][2]
            end

            # Store displacements in variables datax/datay
            datax[cj, ci] -= (max_x1 - M) + idx[cj, ci]
            datay[cj, ci] -= (max_y1 - M) + idy[cj, ci]
            xx[cj, ci] = ii + M / 2
            yy[cj, ci] = ii + N / 2
            ci += 1
        end

        cj += 1
    end
    return xx, yy, datax, datay
end

"""
    finalpass(A::Matrix{Float32}, B::Matrix{Float32}, N::Int32, ol::Float32, 
                    idx::Matrix{Float32}, idy::Matrix{Float32}, Dt::Int32)

    Provides the final pass to get the displacements with
    subpixel resolution.
    
    Parameters
    ----------
    A: Matrix containing image data of first frame.
    B: Matrix containing image data of second frame.
    N: Vector representing sub-window pixel sizes for the pass.
    overlap: Fraction of window overlap. Int.
    idx: Matrix of same type as A, containing data displacement information.
    idy: Matrix of same type as A, containing data displacement information.

    1999 - 2011, J. Kristian Sveen (jks@math.uio.no)
    For use with MatPIV 1.7, Copyright
    Distributed under the terms of the GNU - GPL license
    timestamp: 09:26, 4 Mar 2011
"""
function finalpass(A::Matrix{T}, B::Matrix{T}, N::Int32, ol::Float32,
    idx::Matrix{T}, idy::Matrix{T}, Dt::Int32) where {T}
    M = N

    # FFT setup
    pad_matrix_a = pad_for_xcorr(A[1:M, 1:N])
    pad_matrix_b = pad_for_xcorr(B[1:M, 1:N])
    P = plan_fft(pad_matrix_a; flags=FFTW.ESTIMATE)
    Pi = plan_ifft(pad_matrix_a; flags=FFTW.ESTIMATE)

    # Preallocations
    size_R_1 = size(pad_matrix_a, 1) - 1
    max_x1::Int32 = 0
    max_y1::Int32 = 0
    sy, sx = size(A)
    dim_1 = ceil(Int32, (sy - N) / ((1 - ol) * N)) + 1
    dim_2 = ceil(Int32, (sx - M) / ((1 - ol) * M)) + 1
    up = zeros(eltype(A), (dim_1, dim_2))
    vp = zeros(eltype(A), (dim_1, dim_2))
    SnR = zeros(eltype(A), (dim_1, dim_2))
    Pkh = zeros(eltype(A), (dim_1, dim_2))
    xp = zeros(eltype(A), (dim_1, dim_2))
    yp = zeros(eltype(A), (dim_1, dim_2))

    cj = 1
    # Main pass loop
    for jj in 1:((1-ol)*N):sy-N+1
        ci = 1
        for ii in 1:((1-ol)*M):sx-M+1
            if isnan(idx[cj, ci])
                idx[cj, ci] = 0
            end

            if isnan(idy[cj, ci])
                idy[cj, ci] = 0
            end

            if (jj + idy[cj, ci]) < 1
                idy[cj, ci] = 1 - jj
            elseif (jj + idy[cj, ci]) > (sy - N + 1)
                idy[cj, ci] = sy - N + 1 - jj
            end

            if (ii + idx[cj, ci]) < 1
                idx[cj, ci] = 1 - ii
            elseif (ii + idx[cj, ci]) > (sx - M + 1)
                idx[cj, ci] = sx - M + 1 - ii
            end

            D2 = B[Int32(jj + idy[cj, ci]):Int32(jj + N - 1 + idy[cj, ci]),
                   Int32(ii + idx[cj, ci]):Int32(ii + M - 1 + idx[cj, ci])]
            E = A[Int32(jj):Int32(jj + N - 1), Int32(ii):Int32(ii + M - 1)]

            stad1 = std(E)
            stad2 = std(D2)

            if stad1 == 0
                stad1 = 1
            end
            if stad2 == 0
                stad2 = 1
            end

            E = E .- mean(E)
            F = D2 .- mean(D2)

            # Cross correlate and FFT
            R::Matrix{Float32} = xcorrf2(E, F, P, Pi, pad_matrix_a, pad_matrix_b) ./ (N * M * stad1 * stad2)

            if !any(isnan.(R)) & !all(x -> x == 0, R)
                max_coords = Vector{Tuple{Float32,Float32}}()
                # Find position of maximal value of R
                if size(R, 1) == (N - 1)
                    fast_max!(max_coords, R)
                else
                    subset = R[Int32(0.5 * N + 2):Int32(1.5 * N - 3),
                        Int32(0.5 * M + 2):Int32(1.5 * M - 3)]
                    fast_max!(max_coords, subset)

                    # Adjust indices for subset
                    max_coords = [(i[1] + Int32(0.5 * N + 1), i[2] + Int32(0.5 * M + 1))
                                  for i in max_coords]
                end

                # Handle a vector that has multiple maximum coordinates.
                # Sum the product of each x and y indice with its own indice within
                # the max_coords vector.
                if length(max_coords) > 1
                    max_x1 = round(Int32,
                        sum([c[2]^2 for c in max_coords]) /
                        sum([c[2] for c in max_coords]))
                    max_y1 = round(Int32,
                        sum([c[1]^2 for c in max_coords]) /
                        sum([c[1] for c in max_coords]))
                end

                # Unpack cartesian index type. Only a handful iterations here
                if length(max_coords) == 1
                    max_y1, max_x1 = max_coords[1][1], max_coords[1][2]
                end

                # Some kind of manual adjustment?
                if max_x1 == 1
                    max_x1 = 2
                end
                if max_y1 == 1
                    max_y1 = 2
                end

                # Runs without error, seems to be equivalent.
                # 3-point peak fit using gaussian fit
                x_0, y_0 = intpeak(max_x1, max_y1,
                    R[max_y1, max_x1],
                    R[max_y1, max_x1-1],
                    R[max_y1, max_x1+1],
                    R[max_y1-1, max_x1],
                    R[max_y1+1, max_x1],
                    N
                )

                R2 = copy(R)

                # This section had a note to try to simplify their try-catch
                # clause by using a distance check.
                try
                    R2[max_y1-3:max_y1+3, max_x1-3:max_x1+3] .= NaN
                catch
                    R2[max_y1-1:max_y1+1, max_x1-1:max_x1+1] .= NaN
                end

                # if size(R, 1) == (N - 1)
                if size_R_1 == (N - 1)
                    max_val = maximum(R2)
                    p2_coords = findall(x -> x == max_val, R2)
                    if length(p2_coords) == 1
                        p2_y2, p2_x2 = p2_coords[1][1], p2_coords[1][2]
                    end
                else
                    # Find subset of R2 and remember where all !NaN vals are for later reindexing
                    subset = R2[Int32(0.5 * N):Int32(1.5 * N - 1),
                        Int32(0.5 * M):Int32(1.5 * M - 1)]
                    non_nans = findall(x -> !isnan(x), subset)

                    # Filter out all non-NaN values to find the maximum
                    # NaN vals contaminate arrays for the min/max functions
                    filt_subset = filter(!isnan, subset)
                    subset_max_val = maximum(filt_subset)
                    max_coords_2 = findall(x -> x == subset_max_val, filt_subset)

                    # Re-adjust coordinates to match original coords in R2
                    og_max_coords = non_nans[max_coords_2]
                    max_coords_2 = [(i[1] + floor(Int32, 0.5 * N - 1),
                        i[2] + floor(Int32, 0.5 * M - 1))
                                    for i in og_max_coords]
                end

                if length(max_coords_2) == 1
                    p2_y2, p2_x2 = max_coords_2[1][1], max_coords_2[1][2]
                elseif length(max_coords_2) > 1
                    p2_x2 = round(Int32, length(max_coords_2) ./ 2)
                    p2_y2 = round(Int32, length(max_coords_2) ./ 2)
                elseif isempty(max_coords)
                    error("Empty set found in Final Pass")
                end

                snr = R[max_y1, max_x1] / R2[p2_y2, p2_x2]
                SnR[cj, ci] = snr
                up[cj, ci] = (-x_0 + idx[cj, ci]) / Dt
                vp[cj, ci] = (-y_0 + idy[cj, ci]) / Dt
                Pkh[cj, ci] = R[max_y1, max_x1]
                xp[cj, ci] = ii + (M / 2) - 1
                yp[cj, ci] = ii + (M / 2) - 1

            else
                up[cj, ci] = NaN
                vp[cj, ci] = NaN
                SnR[cj, ci] = NaN
                Pkh[cj, ci] = 0
                xp[cj, ci] = ii + (M / 2) - 1
                yp[cj, ci] = ii + (M / 2) - 1
            end
            ci += 1
        end
        cj += 1
    end
    return xp, yp, up, vp, SnR, Pkh
end


# FOURIER
"""
    pad_for_xcorr(trunc_matrix::Matrix{Float32})

    Determine the dimensions of a padded matrix to be used for xcorrf2. Depends
    on the dimensions of the current window size and the original array. Takes
    the current window size, scales it by a power of 2, then creates an array of 
    zeros that size. 
    
    Parameters
    ----------
        trunc_matrix: The original matrix to be padded, but truncated down to 
            the window's size. 

    Returns
    ----------
        A padded matrix of zeros up to the next power of 2 of the window size.
"""
function pad_for_xcorr(trunc_matrix::Matrix{Float32})
    ma, na = size(trunc_matrix)
    mf = nextpow(2, ma + na)
    return zeros(eltype(ComplexF32), mf, mf)
end


"""
    xcorrf2(A, B, plan, iplan, pad_matrix_a, pad_matrix_b)

    Two-dimensional cross-correlation using Fourier transforms.
    XCORRF2(A,B) computes the crosscorrelation of matrices A and B.

    Parameters
    ----------
    A: matrix (2D array) to be compared.
    B: matrix ((2D array)) to be compared.
    plan: Callable function representing a pre-planned, omptimized FFT. 
        Created using the dimensions of matrices A & B.
    iplan: Callable function representing a pre-planned, optimized inverse FFT.
    
    Returns
    ----------
    c: A matrix whose values reflect the 2D correlation between every cell in 
        A & B.

    Originally written in Matlab by,
    Author(s): R. Johnson
    Revision: 1.0   Date: 1995/11/27
"""
function xcorrf2(A::Matrix{Float32}, B::Matrix{Float32}, 
    plan::FFTW.cFFTWPlan, iplan::AbstractFFTs.ScaledPlan,
    pad_matrix_a::Matrix{ComplexF32}, pad_matrix_b::Matrix{ComplexF32})

    # Unpack size() return tuple into appropriate variables
    ma, na = Int32.(size(A))
    mb, nb = Int32.(size(B))

    # Reverse conjugate
    B_conj = @view B[mb:-1:1, nb:-1:1]
    for i in 1:mb, j in 1:nb
        pad_matrix_b[i, j] = conj(B_conj[i, j])
    end

    # Transfer data from og matrix to optimized sized ones
    pad_matrix_a[1:ma, 1:na] = A[1:ma, 1:na]

    # FFT and trim
    return real(iplan * ((plan * pad_matrix_b) .* (plan * pad_matrix_a))
    )[1:ma+mb-1, 1:na+nb-1]
end


# FILTERS and Interpolations
"""
### linear_naninterp
    Interpolates NaN's in a vectorfield. Sorts all spurious 
    vectors based on the number of spurious neighbors to a 
    point. The function replaces NaN values in the input 
    vector fields U and V using linear interpolation.\n
    Interpolation starts with the ones that have the least 
    number of outliers in their neighborhood and loops until no 
    NaN's are present in the field.\n
    NOTE: This function is a completely gutted version of the
    original naninterp combined with naninterp2. The majority
    of naninterp was varargin code figuring out which method
    to execute naninterp2 with. Since arguments to naninterp
    were hardcoded in, the desired form of naninterp was just
    turned into a specific function. If necessary, other 
    methods can be written to emulate the functionality of
    the original functions. Moreover, though a mask was being
    passed in, the masking branch was not being executed in
    either functions. Therefore, in this implementation, masking
    was left out.\n

    Parameters:
    -----------
    - u, v: `Matrices`

    Original author:
    ----------------
    J. Kristian Sveen (jks@math.uio.no)
    Department of Mathematics, Mechanics Division, University of Oslo, Norway
    Copyright 1999 - 2001
    For use with MatPIV 1.6, Copyright
    Distributed under the terms of the GNU - GPL license
"""
function linear_naninterp(u::Matrix{Float32}, v::Matrix{Float32})
    coords::Vector{CartesianIndex{2}} = findall(x -> isnan(x), u)
    numm = length(coords)
    dy, dx = size(u)
    lp = 1
    tel = 1

    # Now sort the NaN's after how many neighbors they have that are
    # physical values. Then we first interpolate those that have 8
    # neighbors, followed by 7, 6, 5, 4, 3, 2 and 1
    while !isempty(coords)
        nei = zeros(Int32, length(coords), 3)

        # Check neighbors
        for i in eachindex(coords)
            py = coords[i][1]
            px = coords[i][2]
            corx1 = 0
            corx2 = 0
            cory1 = 0
            cory2 = 0

            # Correct if vector is on edge of matrix
            # These are edge cases as we explore each NaN in u.
            if py == 1
                cory1 = 1
                cory2 = 0
            elseif py == dy
                cory1 = 0
                cory2 = -1
            end
            if px == 1
                corx1 = 1
                corx2 = 0
            elseif px == dx
                corx1 = 1
                corx2 = -1
            end

            # Create a matrix of NaN's 8 neighbors
            ma::Matrix{Float32} = u[
                py-1+cory1:py+1+cory2,
                px-1+corx1:px+1+corx2
            ]

            nei[i, 1] = count(!isnan, ma)
            nei[i, 2] = px
            nei[i, 3] = py
        end

        # Sort NEI by row to interpolate vectors with fewest spurious neighbors. -------------------------- CAN COMMENT OUT FOR PERFORMANCE. LOWERS OVERALL QUALITY A BIT THOUGH
        nei::Matrix{Int32} = sortslices(nei, dims=1, lt=Base.isgreater)

        # Reconstruct sorted outliers and interpolate 1st 50%.
        idx = findall(x -> x >= 8, nei[:, 1])
        while isempty(idx)
            idx = findall(x -> x >= (8 - tel), nei[:, 1])
            tel += 1
        end
        tel = 1
        py = nei[idx, 3]
        px = nei[idx, 2]

        for j in axes(py, 1)
            corx1 = 0
            corx2 = 0
            cory1 = 0
            cory2 = 0
            if py[j] == 1
                cory1 = 1
                cory2 = 0
            elseif py[j] == dy
                cory1 = 0
                cory2 = -1
            end
            if px[j] == 1
                corx1 = 1
                corx2 = 0
            elseif px[j] == dx
                corx1 = 0
                corx2 = -1
            end

            # Super hard to read but nice and fancy fast way to take the mean
            # of all the non-NaN values in the area of interest.
            u[py[j], px[j]] = mapreduce(x -> isnan(x) ? 0 : x, +,
                u[py[j]-1+cory1:py[j]+1+cory2,
                    px[j]-1+corx1:px[j]+1+corx2]) /
                              count(!isnan, u[py[j]-1+cory1:py[j]+1+cory2,
                px[j]-1+corx1:px[j]+1+corx2])
            v[py[j], px[j]] = mapreduce(x -> isnan(x) ? 0 : x, +,
                v[py[j]-1+cory1:py[j]+1+cory2,
                    px[j]-1+corx1:px[j]+1+corx2]) /
                              count(!isnan, v[py[j]-1+cory1:py[j]+1+cory2,
                px[j]-1+corx1:px[j]+1+corx2])

            if lp > numm
                u[py[j], px[j]] = 0
                v[py[j], px[j]] = 0
            end
        end
        coords = findall(x -> isnan(x), u)  # Might not be necessary
        lp += 1
    end
    return u, v
end


"""
        regular_interp(samples, xs, ys, XI, YI)

    Interpolates the given `samples` using linear interpolation on a regular grid defined by `xs` and `ys`.
    The interpolated values are evaluated at the points defined by `XI` and `YI`.

    Arguments
    ---------
    - `samples::AbstractArray`: The array of samples to be interpolated.
    - `xs::AbstractArray`: The array of x-coordinates defining the grid.
    - `ys::AbstractArray`: The array of y-coordinates defining the grid.
    - `XI::AbstractArray`: The array of x-coordinates where the interpolated values are evaluated.
    - `YI::AbstractArray`: The array of y-coordinates where the interpolated values are evaluated.

    Returns
    ---------
    - `itp_results`: The interpolated values evaluated at the points defined by `XI` and `YI`.
"""
function regular_interp(samples::Matrix{Float32}, xs::T, ys::T, XI::T, YI::T) where {T}
    try
        itp = Interpolations.interpolate((ys, xs), samples, Gridded(Linear()))
        # Add 0 where extrapolation should be
        extp = Interpolations.extrapolate(itp, 0)  
        # itp_results = [itp(yi, xi) for yi in YI, xi in XI]

        itp_results = zeros(Float32, (size(YI, 1), size(XI, 1)))
        for (j, yi) in enumerate(YI), (i, xi) in enumerate(XI)
            # Extrapolate along the outside border
            if yi == YI[1] || xi == XI[1] || yi == YI[end] || xi == XI[end]
                itp_results[j, i] = extp(yi, xi)
            else
                # Otherwise, justi interpolate as usual
                itp_results[j, i] = itp(yi, xi)
            end
        end
        return itp_results

    catch
        ys = collect(ys)  # Convert to an array to fit into deduplicate_knots!
        xs = collect(xs)
        knots = (ys, xs)
        Interpolations.deduplicate_knots!(knots)
        itp = Interpolations.interpolate(knots, samples, Gridded(Linear()))
        itp_results = [itp(yi, xi) for yi in YI, xi in XI]
        return itp_results
    end
end

function build_grids(ol::Float32, sy::Int32, sx::Int32, win_size::Int32, final_win::Int32)
    windiv2 = win_size ÷ 2
    coarse_step = (1 - ol) * 2 * windiv2
    fine_step = (1 - ol) * windiv2

    if win_size != final_win
        coarse_xs = windiv2 + 1:coarse_step:(sx - windiv2) + 1
        coarse_ys = windiv2 + 1:coarse_step:(sy - windiv2) + 1
        fine_XI = (windiv2 ÷ 2) + 1:fine_step: (sx - windiv2 ÷ 2) + 1 
        fine_YI = (windiv2 ÷ 2) + 1:fine_step: (sy - windiv2 ÷ 2) + 1
    else
        # Will need adjusting below? Neither mlab nor julia enter this branch?
        fine_XI = 1:fine_step: sx - (windiv2 + 1) + windiv2 ÷ 2
        fine_YI = 1:fine_step: sy - (windiv2 + 1) + windiv2 ÷ 2
        coarse_xs = copy(fine_XI)
        coarse_ys = copy(fine_YI)
    end

    # Adjust for odd crops
    if coarse_xs[end] < (sx - windiv2)
        coarse_xs = windiv2 + 1:coarse_step:(sx - windiv2) + windiv2 + 1 
    end
    if coarse_ys[end] < (sy - windiv2)
        coarse_ys = windiv2 + 1:coarse_step:(sy - windiv2) + windiv2 + 1
    end

    return coarse_ys, coarse_xs, fine_YI, fine_XI
end
"""
        build_grids_2(data)

    Builds coarse and fine grids based on the given data.

    Arguments
    ---------
    `data`: A 2-dimensional array representing the data.

    Returns
    ---------
    `coarse_ys`: A 1-dimensional array representing the coarse grid in the y-direction.
    `coarse_xs`: A 1-dimensional array representing the coarse grid in the x-direction.
    `fine_YI`: A 1-dimensional array representing the fine grid in the y-direction.
    `fine_XI`: A 1-dimensional array representing the fine grid in the x-direction.
"""
function build_grids_2(data::Matrix{Float32})

    coarse_y_dim = size(data, 1)
    coarse_x_dim = size(data, 2)

    min_y, max_y = minimum(data[:, 1]), maximum(data[:, 1])
    min_x, max_x = minimum(data[1, :]), maximum(data[1, :])

    # Avoids grid duplication error by adding a slight variation to the data.
    if max_y == min_y
        max_y += Float32(1e-10)
    end
    if max_x == min_x
        max_x += Float32(1e-10)
    end

    coarse_ys = LinRange(min_y, max_y, coarse_y_dim)
    coarse_xs = LinRange(min_x, max_x, coarse_x_dim)

    # NaN border adjustment
    fine_yi_dim = (coarse_y_dim * 2) - 1
    fine_xi_dim = (coarse_x_dim * 2) - 1

    fine_YI = LinRange(min_y, max_y, fine_yi_dim)
    fine_XI = LinRange(min_x, max_x, fine_xi_dim)

    return coarse_ys, coarse_xs, fine_YI, fine_XI
end

"""
    make_border(data::Matrix{Float32})

    Create a bordered matrix around the given data using the outside values of
    the matrix. Copies the 1st/last rows and columns to the new matrix, then
    fills in the interior with the original data.

    Parameters
    ----------
        data: A 2-dimensional array representing the data.

    Returns
    ----------
        bordered_matrix: A 2-dimensional array representing the original matrix
            with a border around it.
"""
function make_border(data::Matrix{Float32})
    # Allocate space for the bordered matrix
    bordered_matrix = zeros(Float64, (size(data, 1) + 2, size(data, 2) + 2))

    # Extract the first and last rows and columns
    first_row = data[1, :]
    last_row = data[end, :]
    first_col = data[:, 1]
    last_col = data[:, end]

    # Place border around new matrix, keeping space for each upcoming new row
    # and column
    bordered_matrix[1, 2:end-1] = first_row
    bordered_matrix[end, 2:end-1] = last_row
    bordered_matrix[2:end-1, 1] = first_col
    bordered_matrix[2:end-1, end] = last_col

    # Fill in the interior of the matrix
    bordered_matrix[2:end-1, 2:end-1] = data
    return bordered_matrix
end


"""
    make_nan_border(data::Matrix{Float32})

    Create a bordered matrix around the given data using NaN values. Fills the
    center of the NaN matrix with the original data.
    
    Parameters
    ----------
        data: A 2-dimensional array representing the data.

    Returns
    ----------
        bordered_matrix: A 2-dimensional array representing the original matrix
            with a NaN border around it one element thick.
"""
function make_nan_border(data::Matrix{Float32})
    # Allocate space for the bordered matrix, make it NaNs then replace them
    bordered_matrix::Matrix{Float32} = fill(NaN, (size(data, 1) + 2, size(data, 2) + 2))

    # Fill in the interior of the matrix
    bordered_matrix[2:end-1, 2:end-1] = data
    return bordered_matrix
end

"""
    globfilt(u::Matrix{Float32}, v::Matrix{Float32})

    Global histogram operator. Find the maximum and minimum velocities allowed 
    vector fields u and v. 

    Arguments
    ---------
        `u`: A 2-dimensional array representing the u-component of the vector field.
        `v`: A 2-dimensional array representing the v-component of the vector field.

    Returns
    ---------
        `u`: A 2-dimensional array representing the filtered u-component.
        `v`: A 2-dimensional array representing the filtered v-component.

    Original Author:
    1999 -2014 copyright J.K.Sveen jks@math.uio.no
    Dept. of Mathematics, Mechanics Division, University of Oslo, Norway
"""
function globfilt(u::Matrix{Float32}, v::Matrix{Float32})
    # println("Global filter running - with limit: 3 * std [U V]")
    # Opportunity for perfomance improvement here.
    nan_filt_u = filter(!isnan, u)
    nan_filt_v = filter(!isnan, v)
    scale_test = maximum(sqrt.(nan_filt_u[:] .^ 2 + nan_filt_v[:] .^ 2))
    if scale_test != 0
        scale = 2 / scale_test
    else
        scale = 0.1
    end

    sx = 3 * nan_std(u[:])
    sy = 3 * nan_std(v[:])
    xo = nan_mean(u[:])
    yo = nan_mean(v[:])
    ii = [xo + sx; xo + sx; xo - sx; xo - sx]
    jj = [yo + sy; yo - sy; yo - sy; yo + sy]

    # Locate points inside chosen area
    poly_points = Point.(ii, jj)  # Define polygon
    points = vec(Point.(u, v))    # Convert u, v to Points for isinside
    in = [isinside(p, poly_points) for p in points]  # Check if points are inside polygon

    u[.!in] .= NaN
    v[.!in] .= NaN

    return u, v
end

"""
    localfilt(u::Matrix{Float32}, v::Matrix{Float32}, 
                        threshold, median_bool=true, m=3)

    Filter out vectors that deviate from the median or the mean of their 
    surrounding neighbors by the factor `threshold` times the standard deviation
    of the neighbors.

    Parameters:
    -----------
    - u, v : `Matrices`
    - threshold : `Int`
        Specifies the point at which a vector has deviated too 
        far from the specified statistical mean or median.
    - median_bool: `Bool`
        If true, specifies that the median should be the turning
        point for the data to be filtered out on. Defaults
        to true. If specified as false, the mean value will be 
        used instead.
        method.
    - m : `Int`
        Defines the number of vectors contributing to the median 
        or mean value of each vector. Defaults to 3, though the
        original implementation mentions that 5 is a good number
        too. Also known as "kernelsize"

    Returns:
    --------
    - hu, hv : `Matrices`
            Successfully filtered matrices. New versions of u
            and v.
"""
function localfilt(x::Matrix{Float32}, y::Matrix{Float32}, u::Matrix{Float32}, 
                    v::Matrix{Float32}, threshold::Int32, median_bool=true, m=3)
    # IN = zeros(eltype(u), size(u))  # Masking not implemented
    dim1 = round(Int32, size(u, 1) + 2 * floor(m / 2))
    dim2 = round(Int32, size(u, 2) + 2 * floor(m / 2))
    nu = zeros(eltype(u), (dim1, dim2)) * NaN
    nv = zeros(eltype(u), (dim1, dim2)) * NaN

    # Transfer over data
    from_cols = round(Int32, floor(m / 2) + 1)
    minus_rows = round(Int32, floor(m / 2))
    nu[from_cols:end-minus_rows, from_cols:end-minus_rows] = u
    nv[from_cols:end-minus_rows, from_cols:end-minus_rows] = v

    # writedlm("../../tests/gpu_tests/nu.csv", nu, ',')
    # writedlm("../../tests/gpu_tests/nv.csv", nv, ',')


    # INx = zeros(eltype(nu), size(nu))
    # Masking not implemented
    # INx[from_cols:end-minus_rows, from_cols:end-minus_rows] = IN

    U2::Matrix{ComplexF32} = nu .+ im .* nv

    ma, na = size(U2)
    histostd = zeros(ComplexF32, size(nu))
    histo = zeros(ComplexF32, size(nu))

    for ii in m-1:1:na-m+2
        for jj in m-1:1:ma-m+2
            # Get a 3x3 submatrix of U2
            m_floor_two = floor(Int32, m / 2)
            tmp = U2[jj-m_floor_two:jj+m_floor_two,
                        ii-m_floor_two:ii+m_floor_two]

            # Assign the center value to NaN
            tmp[ceil(Int32, m / 2), ceil(Int32, m / 2)] = NaN

            # Run the appropriate stat depending on method arg.
            # usum = median_bool ? im_median_magnitude(tmp[:]) : mean(tmp[:])
            histo[jj, ii] = median_bool ? im_median_magnitude(tmp[:]) : mean(tmp[:])
            histostd[jj, ii] = im_std(tmp[:])
        end
    end

    # Locate gridpoints w/higher value than the threshold
    coords = findall(
        (real(U2) .> real(histo) .+ threshold .* real(histostd)) .|
        (imag(U2) .> imag(histo) .+ threshold .* imag(histostd)) .|
        (real(U2) .< real(histo) .- threshold .* real(histostd)) .|
        (imag(U2) .< imag(histo) .- threshold .* imag(histostd))
        )

    # Then "filter" those points out by changing them to NaN!
    for jj in eachindex(coords)
        nu[coords[jj]] = NaN
        nv[coords[jj]] = NaN
    end

    m_ceil_two = ceil(Int32, m / 2)
    m_floor_two = floor(Int32, m / 2)
    hu::Matrix{Float32} = nu[m_ceil_two:end-m_floor_two, m_ceil_two:end-m_floor_two]
    hv::Matrix{Float32} = nv[m_ceil_two:end-m_floor_two, m_ceil_two:end-m_floor_two]

    return hu, hv
end

"""
    intpeak(x1, y1, R, Rxm1, Rxp1, Rym1, Ryp1, N)

    Interpolation uses Gaussian method.
    Interpolates correlation peaks in PIV. 

    Parameters:
    -----------
        x1, y1 : Maximual values in respective directions
        N : Interrogation window size
        Rxm1, Rxp1 : X-max values in matrix R "minus" or "plus" 1.
        Rym1, Ryp1 : Same as above but Y vals.
        R: Matrix resulting from xcorrf2.

    Original Author:
    Time stamp: 12:32, Apr. 14, 2004.
    Copyright 1998-2004, J. Kristian Sveen, 
    jks@math.uio.no/jks36@damtp.cam.ac.uk
    Dept of Mathmatics, University of Oslo/ 
    DAMTP, Univ. of Cambridge, UK
    Distributed under the GNU general public license.
"""
function intpeak(x1::Int32, y1::Int32, R::T, Rxm1::T, Rxp1::T,
    Rym1::T, Ryp1::T, N::Int32) where {T}
    if length(N) == 2
        M = N[1]
        N = N[2]
    else
        M = N
    end
    x01 = x1 + ((log(Complex(Rxm1)) - log(Complex(Rxp1))) /
                ((2 * log(Complex(Rxm1))) - (4 * log(Complex(R))) + (2 * log(Complex(Rxp1)))))

    y01 = y1 + ((log(Complex(Rym1)) - log(Complex(Ryp1))) /
                ((2 * log(Complex(Rym1))) - (4 * log(Complex(R))) + (2 * log(Complex(Ryp1)))))

    x0 = x01 - M
    y0 = y01 - N

    x0 = real(x0)
    y0 = real(y0)

    return x0, y0
end


# STATISTICS
"""
## im_median_magnitude
    Take the median of a collection of complex numbers using the absolute magnitude.
    This great function was created by the Julia Community, specifically:
    @PeterSimmon & @mbauman 07/23/2024
    https://discourse.julialang.org/t/median-of-complex-numbers-different-from-matlab-output/117352/5
"""
function im_median_magnitude(collection::AbstractArray{ComplexF32})
    i = filter(x -> !isnan(x), collection)
    isempty(i) && return NaN
    n = length(i)
    v = partialsort!(i, div(n + 1, 2, RoundDown):div(n + 1, 2, RoundUp); by=x -> (abs2(x), angle(x)))
    return sum(v) / length(v)
end

"""
### im_mean
    Find the mean of the argued collection of complex numbers.
    If the collection is empty, returns NaN.
"""
function im_mean(collection::AbstractArray{ComplexF32})
    if length(collection) < 1
        return NaN
    end
    real_part = mean(real.(collection))
    im_part = mean(imag.(collection))
    return real_part + im_part * im
end

"""
### im_std
    Find the std dev of the argued collection of complex numbers.
    If the collection is empty, returns NaN.
"""
function im_std(collection::AbstractArray{ComplexF32})
    i = filter(x -> !isnan(x), collection)

    if length(i) > 0
        real_part = std(real(i), corrected=false)
        im_part = std(imag(i), corrected=false)
        return real_part + im_part * im
    end

    # All NaN!
    return NaN
end

"""
        nan_std(collection)

    Compute the standard deviation of a collection, excluding any NaN values.

    Arguments
    ---------
        `collection`: The collection of values.

    Returns
    ---------
        The standard deviation of the collection, excluding NaN values. If the 
        collection is empty or contains only NaN values, NaN is returned.
"""
function nan_std(collection::Vector{Float32})
    i = filter(x -> !isnan(x), collection)
    if length(i) > 0
        return std(i, corrected=false)
    end
    return NaN
end

"""
        nan_mean(collection)

    Compute the mean of a collection, excluding any NaN values.

    Arguments
    ---------
    - `collection`: A collection of values.

    Returns
    --------
    - The mean of the collection, excluding NaN values. If the collection is 
    empty or contains only NaN values, NaN is returned.
"""
function nan_mean(collection::Vector{Float32})
    i = filter(x -> !isnan(x), collection)
    if length(i) > 0
        return mean(i)
    end
    return NaN
end

"""
    fast_max!(max_coords::Vector{Tuple{Int32, Int32}}, collection::Matrix{Float32})

    This function finds the maximum value in the given `collection` matrix and 
    updates the `max_coords` vector with the coordinates of the maximum value. 
    If there are multiple maximum values, all their coordinates are added to the 
    `max_coords` vector.

Arguments
---------
- `max_coords::Vector{Tuple{Int32, Int32}}`: A vector to store the coordinates 
                                                of the maximum value(s).
- `collection::Matrix{Float32}`: The matrix in which to find the maximum value(s).

Returns
--------
- `max_coords::Vector{Tuple{Int32, Int32}}`: The updated vector containing the coordinates of the maximum value(s).
"""
function fast_max!(max_coords::Vector{NTuple{2, Float32}}, collection::Matrix{Float32})
    if all(x -> x == 0, collection)
        empty!(max_coords)
        return
    end

    max_val::Float64 = -Inf
    for i in axes(collection, 1)
        for j in axes(collection, 2)
            if collection[i, j] > max_val
                max_val = collection[i, j]
                empty!(max_coords)
                push!(max_coords, (i, j))
            else
                if collection[i, j] == max_val
                    push!(max_coords, (i, j))
                end
            end
        end
    end
    return max_coords
end

function run_test_data()
    fp1 = "/home/server/pi/homes/shindelr/Nearshore-PIV/juliaPIV/data/A044_C001_020824.037047.jpg"
    fp2 = "/home/server/pi/homes/shindelr/Nearshore-PIV/juliaPIV/data/A044_C001_020824.037048.jpg"
    im1::Matrix{Gray{N0f8}} = load(fp1)
    im2::Matrix{Gray{N0f8}} = load(fp2)
    # crops = (24, 2425, 1, 2048)
    crops = (1, 3072, 1, 2048)  # 2023 jpgs were uncropped
    im1 = im1[crops[3]:crops[4], crops[1]:crops[2]]
    im2 = im2[crops[3]:crops[4], crops[1]:crops[2]]

    main((im1, im2), Int32(16), Float32(0.5))
end