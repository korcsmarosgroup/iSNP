apt-get -y install python3.7 python3.7-dev python3-setuptools python3-pip
update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.7 1
pip3 install --upgrade numpy
pip3 install --upgrade scipy
pip3 install --upgrade pandas
pip3 install --upgrade matplotlib
pip3 install --upgrade mlxtend
pip3 install --upgrade xgboost
pip3 install --upgrade scikit-learn
pip3 install --upgrade biopython

# python dependencies for the main executor script
pip3 install --upgrade argparse
pip3 install --upgrade minio
pip3 install --upgrade python-qpid-proton
