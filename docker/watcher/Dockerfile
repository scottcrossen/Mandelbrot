FROM ubuntu:16.04

RUN apt-get update \
  && apt-get install -y inotify-tools make perl gcc net-tools curl python \
  && rm -rf /var/lib/apt/lists/*

COPY docker-entrypoint.sh /docker-entrypoint.sh

VOLUME /app
WORKDIR /app

ENTRYPOINT ["/docker-entrypoint.sh"]
CMD []
