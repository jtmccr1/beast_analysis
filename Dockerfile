FROM openjdk:8
# copy the packaged jar file into our docker image
RUN apt-get update && apt-get install --yes --no-install-recommends \
bc

COPY bin/ ./bin/
COPY lib ./lib/
ENV PATH $PATH:bin/
