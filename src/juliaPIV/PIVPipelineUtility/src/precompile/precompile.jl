using PIVPipelineUtility


@info "Precompiling full functionality now..."
PIVPipelineUtility_src_dir = dirname(pathof(PIVPipelineUtility))  # ~/**/.venv/lib64/python3.12/site-packages/juliaPIV/PIVPipelineUtility/src
out = joinpath(PIVPipelineUtility_src_dir, "precompile/out")
in = joinpath(PIVPipelineUtility_src_dir, "precompile/in/batch/precompile.txt")      # Remember that the .jl file takes a txt file as a batch not the dir
try
    PIVPipelineUtility.io_main_wrapper(
        Int32(2),                                           # N
        Int32(24), Int32(2425), Int32(1), Int32(2048),      # Crop Factors
        Int32(16),                                          # Final Window Size
        Float32(0.5),                                       # OL
        out,                                                # Out
        in,                                                 # In
        Int(0),                                             # Quiet?
        Float32(1.0),                                       # Downsample
        Int(0)                                              # Save Images?
    )
catch e
    @error "Error in io_main: $e"
    Base.showerror(stderr, e)
    println(stderr, "\n", sprint(Base.show_backtrace, backtrace()))
end
