"""
Utility to compile the PIVPipelineUtility.jl and main.jl scripts.
"""

import subprocess
import os
import pathlib

def compile_juliapiv():
    """
    Compile juliaPIV on the target OS.
    """
    root_dir = pathlib.Path(__file__).resolve().parent   # juliaPIV/src/juliaPIV/
    pivpipelineutility_path = os.path.join(root_dir, "PIVPipelineUtility")
    precompile_dir_path = os.path.join(pivpipelineutility_path, "src/precompile")
    precompile_script_path = os.path.join(precompile_dir_path, "precompile.jl")

    cmmd = [
        'julia',
        f'--project={pivpipelineutility_path}',
        '-e',
        f"""
        using Pkg;
        using PackageCompiler;
        Pkg.instantiate();
        create_library("{pivpipelineutility_path}", "src/juliaPIV/piv_build"; lib_name="pivbuild", precompile_execution_file="{precompile_script_path}", force=true);
        """
        # isdir("src/juliaPIV/piv_build/lib/julia") || mkdir("src/juliaPIV/piv_build/lib/julia");
        # run(`install_name_tool -add_rpath @loader_path/julia src/juliaPIV/piv_build/lib/libpivbuild.dylib`);
        # cp(joinpath(Sys.BINDIR, "..", "lib", "libjulia.dylib"), "src/juliaPIV/piv_build/lib/julia/"; force=true);
        # cp(joinpath(Sys.BINDIR, "..", "lib", "libjulia-internal.dylib"), "src/juliaPIV/piv_build/lib/julia/"; force=true);
    ]
    subprocess.run(cmmd)

    print(f"Build location: {os.path.join(root_dir, 'piv_build')}")