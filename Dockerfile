FROM rocker/geospatial:4.2.1

# Install GRASS
# https://bugs.launchpad.net/ubuntu/+source/grass/+bug/1995244
RUN --mount=type=cache,sharing=locked,target=/var/cache/apt \
    --mount=type=cache,sharing=locked,target=/var/lib/apt/lists \
    rm -f /etc/apt/apt.conf.d/docker-clean; \
    echo 'Binary::apt::APT::Keep-Downloaded-Packages "true";' \
        > /etc/apt/apt.conf.d/keep-cache && \
    apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -qy \
        grass-core \
        python3-six

# Install package
WORKDIR /home/rstudio
COPY DESCRIPTION .
RUN Rscript -e 'devtools::install(dependencies=TRUE)'

# Copy files
COPY --chown=rstudio:rstudio . .
