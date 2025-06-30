using PIVPipelineUtility

# Dev code
Base.ARGS = ["2", 
"24,2425,1,2048", 
"16", 
"0.5",
"/Users/robinshindelman/repos/Nearshore-Research/tests/abs-batch-testing/piv_out/",
"/Users/robinshindelman/repos/Nearshore-Research/tests/abs-batch-testing/batches/batchtest.txt",
 "1"]

PIVPipelineUtility.julia_main()

println("Precompile complete")
