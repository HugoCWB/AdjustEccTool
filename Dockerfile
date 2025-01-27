# This Dockerfile uses neuropythy's docker image as a starting point
# because the code contained in it uses neuropythy.

# TODO: Change this to match the version tag before publication.
# (Later: make it 0.12.16).
FROM nben/neuropythy:latest

# Note the Maintainer. TODO: add email address.
MAINTAINER Hugo T. Chow-Wing-Bom <hugo.chow-wing-bom.15@ucl.ac.uk>


# Copy our library into the docker image.
USER root
COPY ./ $HOME/repo/
RUN chown -R $NB_USER $HOME/repo


# Install the library locally.
USER $NB_USER
RUN cd $HOME/repo \
 && eval "$(command conda shell.bash hook)" \
 && conda activate \
 && pip install -e .


# Add the entrypoint.
ENTRYPOINT ["python", "-m", "eccen_adjust"]
