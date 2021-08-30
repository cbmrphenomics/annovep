TAG := v1

all:
	docker build -t annovep:${TAG} .
	docker build -t annovep:latest .
	docker image prune --force
