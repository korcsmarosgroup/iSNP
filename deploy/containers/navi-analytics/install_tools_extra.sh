#!/bin/bash
set -e

# install python environment
add-apt-repository ppa:deadsnakes/ppa
apt-get update

# install bioinformatics tools
apt-get -y install curl libz-dev autoconf autogen pkg-config
mkdir /build
cd /build

# RSAT
wget https://korcsmaroslab.org/rsat_installation/rsat_files.tar.gz
tar -xzvf rsat_files.tar.gz
cp /build/rsat_modules/compare-matrices-quick /usr/bin/compare-matrices-quick
cp /build/rsat_modules/count-words /usr/bin/count-words
cp /build/rsat_modules/info-gibbs /usr/bin/info-gibbs
cp /build/rsat_modules/matrix-scan-quick /usr/bin/matrix-scan-quick
cp /build/rsat_modules/retrieve-variation-seq /usr/bin/retrieve-variation-seq
mkdir /rsat
mkdir /rsat/perl-scripts
mkdir /rsat/build
cp -r /build/rsat_perl-scripts/rsat/perl-scripts/* /rsat/perl-scripts

# BEDTOOLS
curl -L https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz | tar xvz
( cd bedtools2 && make )
cp bedtools2/bin/bedtools /usr/bin/bedtools

# VCFTOOLS
wget https://github.com/vcftools/vcftools/archive/v0.1.16.tar.gz
tar -xzvf v0.1.16.tar.gz
( cd vcftools-0.1.16 && ./autogen.sh && ./configure && make )
cp vcftools-0.1.16/src/cpp/vcftools /usr/bin/vcftools

# SNPEFF
wget https://korcsmaroslab.org/snpeff_installation/snpEff.jar
cp snpEff.jar /usr/bin/snpEff.jar
echo "#!/bin/bash\njava -jar snpEff.jar\n" > /usr/bin/snpEff
chmod +x /usr/bin/snpEff
wget https://korcsmaroslab.org/snpeff_installation/SnpSift.jar
cp SnpSift.jar /usr/bin/SnpSift.jar
echo "#!/bin/bash\njava -jar SnpSift.jar\n" > /usr/bin/SnpSift
chmod +x /usr/bin/SnpSift

# MEME
curl -L http://meme-suite.org/meme-software/5.1.1/meme-5.1.1.tar.gz --output meme-5.1.1.tar.gz
tar zxf meme-5.1.1.tar.gz
cd meme-5.1.1
./configure --prefix=/meme --with-url=http://meme-suite.org/ --enable-build-libxml2 --enable-build-libxslt
make
make install
cd /build

# VIENNARNA
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.14.tar.gz
tar -xzvf ViennaRNA-2.4.14.tar.gz
cd ViennaRNA-2.4.14/
./configure --without-perl --without-python --without-python3
make
export PATH=$PATH:$PWD/src/bin
cd /build

# MIRNADA
# wget https://korcsmaroslab.org/miranda_installation/miranda-master.zip
# unzip miranda-master.zip
# cd miranda-master
# ./configure --host=x86_64-Linux && make && make install
mv /usr/local/bin/miranda.sh /usr/local/bin/miranda

cd /
#rm -rf /build
