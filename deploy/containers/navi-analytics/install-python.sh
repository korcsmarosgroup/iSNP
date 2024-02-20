sudo apt-get update
sudo apt-get install build-essential zlib1g-dev libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libreadline-dev libffi-dev libsqlite3-dev wget libbz2-dev
cd /tmp
wget https://www.python.org/ftp/python/3.7.9/Python-3.7.9.tgz
tar -xf Python-3.7.9.tgz
cd Python-3.7.9
./configure --enable-optimizations
make -j4  # Adjust the number according to the number of CPU cores
apt-get -y install python3-pip

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
