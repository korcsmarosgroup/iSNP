apt update
apt-get install -y build-essential zlib1g-dev libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libreadline-dev libffi-dev libsqlite3-dev libbz2-dev liblzma-dev tk-dev libdb-dev libexpat1-dev liblzma-dev libffi-dev libgdbm-dev libglib2.0-dev libssl-dev
wget https://www.python.org/ftp/python/3.7.4/Python-3.7.4.tgz
tar xzf Python-3.7.4.tgz
cd Python-3.7.4
./configure
make
make install

apt-get -y install python3-pip

pip3 install --upgrade numpy
pip3 install --upgrade scipy
pip3 install --upgrade pandas
pip3 install --upgrade matplotlib
pip3 install --upgrade mlxtend
pip3 install --upgrade xgboost
pip3 install --upgrade scikit-learn
pip3 install --upgrade biopython
pip3 install --upgrade docker

# python dependencies for the main executor script
pip3 install --upgrade argparse
pip3 install --upgrade minio
pip3 install --upgrade python-qpid-proton
