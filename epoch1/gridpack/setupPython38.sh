lcgviews=/cvmfs/sft.cern.ch/lcg/views
lcgreleases=/cvmfs/sft.cern.ch/lcg/releases

# Choose the LCG view

python3view=LCG_99/x86_64-centos7-gcc8-opt
###echo Setting up python3 from ${lcgviews}/${python3view}

# Determine the python release for that LCG view

python3rel=$(\ls -l ${lcgviews}/${python3view}/bin/python3 | awk '{print $NF}')
python3rel=${python3rel%/bin/python3}
###echo Setting up ${python3rel}
python3tag=${python3rel#${lcgreleases}/}
echo Setting up ${python3tag}

# Set all the environment variables that appear in "env | grep LCG"
# after running ". ${lcgviews}/${python3view}/setup.sh"

export PATH=${python3rel}/bin:${PATH}
export LD_LIBRARY_PATH=${python3rel}/lib:${LD_LIBRARY_PATH}
export PYTHONPATH=${python3rel}/lib:${python3rel}/lib/python3.8/site-packages:${PYTHONPATH}
export MANPATH=${python3rel}/share/man:${MANPATH}
export C_INCLUDE_PATH=${python3rel}/include:${C_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${python3rel}/include:${CPLUS_INCLUDE_PATH}

# Determine the six release for that LCG view

sixrel=$(\ls -l ${lcgviews}/${python3view}/lib/python3.8/site-packages/six.py | awk '{print $NF}')
sixrel=${sixrel%/lib/python3.8/site-packages/six.py}
echo Setting up ${sixrel}
export PYTHONPATH=${sixrel}/lib/python3.8/site-packages:${PYTHONPATH}


