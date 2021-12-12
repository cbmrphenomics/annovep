FROM docker.io/ensemblorg/ensembl-vep:release_105.0

# Root needed to install packages; the entrypoint switches to a regular user
USER root

# Install all dependencies
COPY ./scripts/install_dependencies.sh /opt/
RUN bash /opt/install_dependencies.sh

# Install VEP plugins
RUN ./INSTALL.pl --AUTO p --PLUGINS all --CACHEDIR /opt/vep-plugins/

# Install VEP plugin LOFTEE (https://github.com/konradjk/loftee)
COPY ./binaries/loftee-v1.0.3.tar.gz /opt/vep-plugins
RUN cd /opt/vep-plugins && \
    tar xvzf ./loftee-v1.0.3.tar.gz && \
    cp -a loftee-1.0.3/* Plugins/

RUN apt-get update && apt-get install python3-ruamel.yaml
RUN pip3 install aush --no-cache

# Create folder for mounting the (shared) cache
RUN mkdir -p /data/cache && touch /data/cache/not_mounted
# Create folder for user data (i.e. the user's current working directory)
RUN mkdir -p /data/user && touch /data/user/not_mounted

COPY ./pipeline/ /opt/annovep/pipeline/
COPY ./scripts/ /opt/annovep/scripts/

# Mountpoint for the current working directory
WORKDIR /data/user

ENTRYPOINT [ "/usr/bin/python3", "/opt/annovep/scripts/entrypoint.py" ]
