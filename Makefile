TAG := v1

all:
	podman build -t annovep:${TAG} .
	podman build -t annovep:latest .

	# FIXME: Needs podman v3.3.1: https://github.com/containers/podman/issues/10832
	# docker image prune --force
