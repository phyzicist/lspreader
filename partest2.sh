set -x
python partest2a.py
mpiexec -n 4 python partest2b.py
