allp = c("boltons==23.0.0", "brotlipy==0.7.0", "certifi==2023.7.22", 
"cffi==1.15.1", 
"cryptography==41.0.3", "idna==3.4", "jsonpatch==1.32", "jsonpointer==2.1", 
"pluggy==1.0.0", 
"pycosat==0.6.6", "pycparser==2.21", "pyOpenSSL==23.2.0", 
"ruamel.yaml==0.17.21", "ruamel.yaml.clib==0.2.6", 
"zstandard==0.19.0", "adjustText==0.7.3", "anndata==0.9.2", "attrs==23.1.0", 
"Brotli==1.1.0", "cached-property==1.5.2", "certifi==2023.7.22", 
"charset-normalizer==3.3.2", "colorama==0.4.6", "contourpy==1.1.1", 
"cycler==0.12.1", "filelock==3.13.1", "fonttools==4.44.0", "fsspec==2023.10.0", 
"gmpy2==2.1.2", "h5py==3.10.0", "importlib-metadata==6.8.0", 
"importlib-resources==6.1.1", "Jinja2==3.1.2", "joblib==1.3.2", 
"kiwisolver==1.4.5", "kneed==0.8.5", "llvmlite==0.41.1", "MarkupSafe==2.1.3", 
"matplotlib==3.7.3", "mpmath==1.3.0", "munkres==1.1.4", "natsort==8.4.0", 
"networkx==3.1", "numba==0.58.1", "numexpr==2.8.4", "numpy==1.24.4", 
"packaging==23.2", "pandas==2.0.3", "patsy==0.5.3", "Pillow==10.0.1", 
"pip==23.3.1", "platformdirs==4.0.0", "pooch==1.8.0", "py-cpuinfo==9.0.0", 
"pybedtools==0.9.1", "pynndescent==0.5.10", "pyparsing==3.1.1", 
"pysam==0.22.0", "PySocks==1.7.1", "python-dateutil==2.8.2", 
"pytz==2023.3.post1", "requests==2.31.0", "scanpy==1.9.5", "scikit-learn==1.3.2", 
"scikit-misc==0.1.4", "scipy==1.10.1", "seaborn==0.13.0", "session-info==1.0.0", 
"setuptools==68.2.2", "simba==1.2", "six==1.16.0", "statsmodels==0.14.0", 
"stdlib-list==0.9.0", "sympy==1.12", "tables==3.8.0", "tbb==2021.10.0", 
"threadpoolctl==3.2.0", "torch==2.1.0", "torchbiggraph==1.0.1.dev0", 
"tqdm==4.66.1", "typing_extensions==4.8.0", "tzdata==2023.3", 
"umap-learn==0.5.4", "unicodedata2==15.1.0", "urllib3==2.0.7", 
"wheel==0.41.3", "zipp==3.17.0")


# necessary for python module control
bsklenv <- basilisk::BasiliskEnvironment(channels=c("conda-forge", "bioconda"),
  envname = "bsklenv", packages = "simba==1.2", pip = setdiff(allp, "simba==1.2"),
  pkgname = "BiocSIMBA")


#numpy>=1.17.0
#pandas>=1.0,!=1.1 #  required by Anndata
#anndata>=0.7.4
## h5py<3.0.0 #  avoid byte strings but caused building errors
## h5py>=3.4
#scikit-learn>=0.19
#scipy>=1.4
#kneed>=0.7
#seaborn>=0.11
#matplotlib>=3.3
#scikit-misc>=0.1.3
#adjusttext>=0.7.3
#umap-learn>=0.3.0
##plotly>=4.14.0
#pybedtools>=0.8.0
## bedtools>=2.29.0 # not available in pip
#tables


