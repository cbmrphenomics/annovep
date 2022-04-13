TAG := v1

all:
	echo "VERSION = '$(shell git log -1 --format="%h (%ai)")'" > pipeline/_version.py
	podman build -t annovep:${TAG} .

# Needs podman v3.3.1: https://github.com/containers/podman/issues/10832
	podman image prune --force
	buildah rm --all
