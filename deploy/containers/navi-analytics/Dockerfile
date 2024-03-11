
FROM navi-base:1.0.0
MAINTAINER NavigOmix project "szalay.beko.mate@gmail.com"

ADD miranda.sh /usr/local/bin/

ADD install_tools_extra.sh /
RUN /bin/bash /install_tools_extra.sh

ADD install-python.sh /
RUN /bin/bash /install-python.sh

ENV PATH /rsat/perl-scripts:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV RSAT /rsat

ENV PATH /meme/bin:/meme/libexec/meme-5.1.1:$PATH

ADD install-perl-libs.sh /
RUN /bin/bash /install-perl-libs.sh

ADD analytic-modules /analytic-modules
ENV PYTHONPATH /analytic-modules
