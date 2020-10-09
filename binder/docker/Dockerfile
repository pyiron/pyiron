FROM alpine:3.8
MAINTAINER Jan Janssen <janssen@mpie.de>

# Inspired by :
# * https://github.com/jupyter/docker-stacks
# * https://github.com/CognitiveScale/alpine-miniconda
# * https://github.com/show0k/alpine-jupyter-docker
# * https://github.com/datarevenue-berlin/alpine-miniconda

# Install glibc and useful packages
RUN echo "@testing http://nl.alpinelinux.org/alpine/edge/testing" >> /etc/apk/repositories \
    && apk --update add \
    bash \
    curl \
    ca-certificates \
    libstdc++ \
    glib \
    git \
    tini@testing \
    && curl "https://alpine-pkgs.sgerrand.com/sgerrand.rsa.pub" -o /etc/apk/keys/sgerrand.rsa.pub \
    && curl -L "https://github.com/sgerrand/alpine-pkg-glibc/releases/download/2.23-r3/glibc-2.23-r3.apk" -o glibc.apk \
    && apk add glibc.apk \
    && curl -L "https://github.com/sgerrand/alpine-pkg-glibc/releases/download/2.23-r3/glibc-bin-2.23-r3.apk" -o glibc-bin.apk \
    && apk add glibc-bin.apk \
    && curl -L "https://github.com/andyshinn/alpine-pkg-glibc/releases/download/2.25-r0/glibc-i18n-2.25-r0.apk" -o glibc-i18n.apk \
    && apk add --allow-untrusted glibc-i18n.apk \
    && /usr/glibc-compat/bin/localedef -i en_US -f UTF-8 en_US.UTF-8 \
    && /usr/glibc-compat/sbin/ldconfig /lib /usr/glibc/usr/lib \
    && rm -rf glibc*apk /var/cache/apk/*

# Configure environment
ENV CONDA_DIR /srv/conda/envs/notebook/
ENV PATH ${CONDA_DIR}/bin:${PATH}
ENV SHELL /bin/bash
ENV PYIRON_USER pyiron
ENV PYIRON_UID 1000
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8
ENV HOME /home/${PYIRON_USER}

# Openmpi fixes
ENV OMPI_MCA_plm isolated
ENV OMPI_MCA_rmaps_base_oversubscribe yes
ENV OMPI_MCA_btl_vader_single_copy_mechanism none

# Configure Miniconda
ENV MINICONDA_VER 4.8.3
ENV MINICONDA Miniconda3-py38_${MINICONDA_VER}-Linux-x86_64.sh
ENV MINICONDA_URL https://repo.continuum.io/miniconda/${MINICONDA}
ENV MINICONDA_MD5_SUM d63adf39f2c220950a063e0529d4ff74

COPY . ${HOME}

# Create user with UID=1000 and in the 'users' group
RUN adduser -s ${SHELL} --disabled-password --gecos "Default user" -u ${PYIRON_UID} -D ${PYIRON_USER} \
    && mkdir -p ${CONDA_DIR} \
    && chown ${PYIRON_USER} ${CONDA_DIR} \
    && chown -R ${PYIRON_USER} ${HOME}

USER ${PYIRON_USER}

# Install 
RUN cd /tmp \
    && mkdir -p ${CONDA_DIR} \
    && curl -L ${MINICONDA_URL}  -o miniconda.sh \
    && echo "${MINICONDA_MD5_SUM}  miniconda.sh" | md5sum -c - \
    && ${SHELL} miniconda.sh -f -b -p ${CONDA_DIR} \
    && rm miniconda.sh \
    && ${CONDA_DIR}/bin/conda install --yes conda==${MINICONDA_VER} \
    && source ${CONDA_DIR}/bin/activate \
    && conda env update -n root -f "${HOME}/binder/environment.yml" \
    && conda install --yes -c conda-forge nodejs \
    && ${SHELL} ${HOME}/binder/postBuild \
    && ${CONDA_DIR}/bin/pip install --force-reinstall --no-deps --pre pyiron \
    && conda clean --all -y \
    && rm ${HOME}/Dockerfile ${HOME}/deploy.sh

USER root

# Configure container startup as root
WORKDIR ${HOME}/
ENTRYPOINT ["/sbin/tini", "--"]
CMD ["jupyter", "notebook", "--ip", "0.0.0.0"]

# Switch back to pyiron to avoid accidental container runs as root
USER ${PYIRON_USER}
