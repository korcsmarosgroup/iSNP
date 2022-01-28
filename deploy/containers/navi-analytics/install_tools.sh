#!/bin/bash
set -e

# install python environment
add-apt-repository -y ppa:jonathonf/python-3.6
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
curl -L https://github.com/vcftools/vcftools/archive/v0.1.16.tar.gz | tar xvz
( cd vcftools-0.1.16 && ./autogen.sh && ./configure && make )
cp vcftools-0.1.16/src/cpp/vcftools /usr/bin/vcftools

# SNPEFF
curl -L https://kent.dl.sourceforge.net/project/snpeff/snpEff_v4_3t_core.zip --output snpEff_v4_3t_core.zip
unzip snpEff_v4_3t_core.zip
cp snpEff/snpEff.jar /usr/bin/snpEff.jar
echo "#!/bin/bash\njava -jar snpEff.jar\n" > /usr/bin/snpEff
chmod +x /usr/bin/snpEff
cp snpEff/SnpSift.jar /usr/bin/SnpSift.jar
echo "#!/bin/bash\njava -jar SnpSift.jar\n" > /usr/bin/SnpSift
chmod +x /usr/bin/SnpSift

# MIRNADA
curl -L http://cbio.mskcc.org/microrna_data/miRanda-aug2010.tar.gz | tar xzv
cd miRanda-3.3a/
./configure && make && make install

cd /
rm -rf /build
