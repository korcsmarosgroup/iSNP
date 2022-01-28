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
curl http://pedagogix-tagc.univ-mrs.fr/download_rsat/previous_versions/rsat_2018-08-01.tar.gz | tar xvz
( cd rsat/contrib/compare-matrices-quick && make clean && make )
( cd rsat/contrib/count-words && make clean && make )
( cd rsat/contrib/info-gibbs && make clean && make )
( cd rsat/contrib/matrix-scan-quick && make clean && make )
( cd rsat/contrib/retrieve-variation-seq && make clean && make )
cp /build/rsat/contrib/compare-matrices-quick/compare-matrices-quick /usr/bin/compare-matrices-quick
cp /build/rsat/contrib/count-words/count-words /usr/bin/count-words
cp /build/rsat/contrib/info-gibbs/info-gibbs /usr/bin/info-gibbs
cp /build/rsat/contrib/matrix-scan-quick/matrix-scan-quick /usr/bin/matrix-scan-quick
cp /build/rsat/contrib/retrieve-variation-seq/retrieve-variation-seq /usr/bin/retrieve-variation-seq
mkdir /rsat
mkdir /rsat/perl-scripts
mkdir /rsat/build
cp -r /build/rsat/perl-scripts/* /rsat/perl-scripts
#cp -r /build/rsat/build/* /rsat/build

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
curl -L -k https://kent.dl.sourceforge.net/project/snpeff/snpEff_v4_3t_core.zip --output snpEff_v4_3t_core.zip
unzip snpEff_v4_3t_core.zip
cp snpEff/snpEff.jar /usr/bin/snpEff.jar
echo "#!/bin/bash\njava -jar snpEff.jar\n" > /usr/bin/snpEff
chmod +x /usr/bin/snpEff
cp snpEff/SnpSift.jar /usr/bin/SnpSift.jar
echo "#!/bin/bash\njava -jar SnpSift.jar\n" > /usr/bin/SnpSift
chmod +x /usr/bin/SnpSift

# MEME
curl -L http://meme-suite.org/meme-software/5.1.1/meme-5.1.1.tar.gz --output meme-5.1.1.tar.gz
tar zxf meme-5.1.1.tar.gz
cd meme-5.1.1
./configure --prefix=/meme --with-url=http://meme-suite.org/ --enable-build-libxml2 --enable-build-libxslt
make
make install
# export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.1.1:$PATH

# VIENNARNA
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.14.tar.gz
tar -xzvf ViennaRNA-2.4.14.tar.gz
cd ViennaRNA-2.4.14/
./configure --without-perl --without-python --without-python3
make
export PATH=$PATH:$PWD/src/bin
cpan Bio::TreeIO
cpan Statistics::Lite

# MIRNADA
curl -L http://cbio.mskcc.org/microrna_data/miRanda-aug2010.tar.gz | tar xzv
cd miRanda-3.3a/
./configure && make && make install

cd /
rm -rf /build
