using PIVPipelineUtility

println("Precompiling PIV C functionality.")
ptr = PIVPipelineUtility.c_io_main_ptr()
@assert ptr != C_NULL "C pointer precompile failed."

# Dev code
println("Precompiling CLI functionality now.")
Base.ARGS = ["2", 
            "24,2425,1,2048", 
            "16", 
            "0.5",
            "src/PIVPipelineUtility/precompile/in/batch",
            "src/PIVPipelineUtility/precompile/out",
            "0",
            "1.0",
            "false"]
PIVPipelineUtility.julia_main()

println("Full Precompile complete")
