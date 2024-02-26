# using a specific 16.04 version (to avoid frequent updates)
FROM ubuntu:16.04
MAINTAINER NavigOmix project "szalay.beko.mate@gmail.com"

ADD install_base_layer.sh /
ADD nav_dev_key.* /

RUN /bin/bash /install_base_layer.sh

ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

EXPOSE 22

WORKDIR /home/navigator

CMD /usr/sbin/sshd -D -o UseDNS=no -o UsePAM=no

