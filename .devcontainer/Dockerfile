FROM mcr.microsoft.com/vscode/devcontainers/anaconda:dev-3

# Copy environment.yml (if found) to a temp locaition so we update the environment. Also
# copy "noop.txt" so the COPY instruction does not fail if no environment.yml exists.
COPY .ci_support/environment.yml /tmp/conda-tmp/
RUN if [ -f "/tmp/conda-tmp/environment.yml" ]; then  sed 's/=/>=/g' /tmp/conda-tmp/environment.yml > /tmp/conda-tmp/environmentflex.yml; /opt/conda/bin/conda env update -n base -f /tmp/conda-tmp/environmentflex.yml; fi \
    && rm -rf /tmp/conda-tmp

# [Optional] Uncomment this section to install additional OS packages.
# RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
#     && apt-get -y install --no-install-recommends <your-package-list-here>
