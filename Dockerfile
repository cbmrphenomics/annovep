FROM ensemblorg/ensembl-vep:release_104.3

# Root needed to install packages; the entrypoint switches to a regular user
USER root

# Install all dependencies
COPY src/install_dependencies.sh /opt/
RUN bash /opt/install_dependencies.sh

# Install VEP plugins
RUN ./INSTALL.pl --AUTO p --PLUGINS all --CACHEDIR /opt/vep-plugins/

# Install VEP plugin LOFTEE (https://github.com/konradjk/loftee)
COPY ./binaries/loftee-v1.0.3.tar.gz /opt/vep-plugins
RUN cd /opt/vep-plugins && \
    tar xvzf ./loftee-v1.0.3.tar.gz && \
    cp -a loftee-1.0.3/* Plugins/

# Install annovar (https://annovar.openbioinformatics.org/)
COPY ./binaries/annovar.scripts.20200608.tar.gz /opt/annovar/
RUN cd /opt/annovar/ && tar xvzf annovar.scripts.20200608.tar.gz --strip-component=1

# Apply tweaks to scripts
COPY ./binaries/annovar.scripts.20200608.patch /opt/annovar/annovar.patch
RUN cd /opt/annovar/ && patch -p1 <annovar.patch

# Create folder for mounting the (shared) cache
RUN mkdir -p /data/cache && touch /data/cache/not_mounted
# Create folder for user data (i.e. the user's current working directory)
RUN mkdir -p /data/user && touch /data/user/not_mounted

COPY ./src/ /opt/annovep/

WORKDIR /data/user

ENTRYPOINT [ "/usr/bin/python3", "/opt/annovep/entrypoint.py" ]
