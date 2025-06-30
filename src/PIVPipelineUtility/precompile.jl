using PIVPipelineUtility

# Dev code
Base.ARGS = ["2", 
"24,2425,1,2048", 
"16", 
"0.5",
# "../tests/pipeline_utility_testing/SVSout_23227179_1724441851/precompile_pivframes",
# "../tests/pipeline_utility_testing/testbatches/tmp.2YNbHiPCwK.txt", 
"/Users/robinshindelman/repos/Nearshore-Research/tests/abs-batch-testing/piv_out/",
"/Users/robinshindelman/repos/Nearshore-Research/tests/abs-batch-testing/batches/batchtest.txt",
 "1"]

PIVPipelineUtility.julia_main()

println("Precompile complete")
