.PHONY:

# Required System files
DOCKER_COMPOSE_EXE := $(shell which docker-compose)

# Variables
ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
MY_UID := $$(id -u)
MY_GID := $$(id -g)
THIS_USER := $(MY_UID):$(MY_GID)

# STDOUT Formatting
RED := $$(echo  "\033[0;31m")
YELLOW := $$(echo "\033[0;33m")
END := $$(echo  "\033[0m")
ERROR_HEADER :=  [ERROR]:
INFO_HEADER := "**************** "
DONE_MESSAGE := $(YELLOW)$(INFO_HEADER) "- done\n" $(END)

# Paths
SCRATCH_DIR := $(ROOT_DIR)/scratch/

# Commands
DOCKER_COMPOSE_CMD := MY_UID=$(MY_UID) MY_GID=$(MY_GID) $(DOCKER_COMPOSE_EXE) -f $(ROOT_DIR)/docker-compose.yml

#############################################################
#  Cleaning targets
#############################################################


#############################################################
#  Building targets
#############################################################

install-dependencies-ubuntu:
	@sudo apt-get install docker.io curl
	@sudo curl -L "https://github.com/docker/compose/releases/download/1.25.4/docker-compose-$$(uname -s)-$$(uname -m)" -o /usr/local/bin/docker-compose
	@sudo chmod +x /usr/local/bin/docker-compose
	@sudo groupadd docker
	@sudo gpasswd -a $$USER docker

#############################################################
#  Docker targets
#############################################################

run-execution:
	@$(DC_UP_CMD) guidance

