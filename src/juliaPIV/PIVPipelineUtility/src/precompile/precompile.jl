using PIVPipelineUtility

@info "Precompiling PIV C functionality..."
ptr = PIVPipelineUtility.c_io_main_ptr()
@assert ptr != C_NULL "C pointer precompile failed."

@info "Precompiling full functionality now..."
out = joinpath(@__DIR__, "out")
in = joinpath(@__DIR__, "in/batch/precompile.txt")      # Remember that the .jl file takes a txt file as a batch not the dir
cstr_out = Base.unsafe_convert(Cstring, out)
cstr_in = Base.unsafe_convert(Cstring, in)
try
    PIVPipelineUtility.c_io_main(
        Int32(2),                                           # N
        Int32(24), Int32(2425), Int32(1), Int32(2048),      # Crop Factors
        Int32(16),                                          # Final Window Size
        Float32(0.5),                                       # OL
        cstr_out,                                           # Out
        cstr_in,                                            # In
        Cint(0),                                            # Quiet?
        Float32(1.0),                                       # Downsample
        Cint(0)                                             # Save Images?
    )
catch e
    @error "Error in io_main: $e"
    Base.showerror(stderr, e)
    println(stderr, "\n", sprint(Base.show_backtrace, backtrace()))
end
