using Pkg
Pkg.add("Conda")
ENV["PYTHON"] = ""
Pkg.add("PyCall")
Pkg.build("PyCall")
using Conda
Conda.add(["batman-package"]; channel="conda-forge")
