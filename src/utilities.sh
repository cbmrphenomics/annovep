#!/bin/bash

if [ -t 1 ]; then
	# STDOUT is a TTY
	readonly TTY_GREEN="\\x1b[32m"
	readonly TTY_GREY="\\x1b[1;30m"
	readonly TTY_RED="\\x1b[31m"
	readonly TTY_YELLOW="\\x1b[33m"
	readonly TTY_END="\x1b[0m"
else
	readonly TTY_GREEN=""
	readonly TTY_GREY=""
	readonly TTY_RED=""
	readonly TTY_YELLOW=""
	readonly TTY_END=""
fi

function _color_log() {
	local -r level="${1}"
	shift 1

	echo -e "${TTY_GREEN}$(date +"%Y-%m-%d %H:%M:%S")${TTY_END} ${TTY_GREY}${level}${TTY_END} $@"
}

function info() {
	_color_log "INFO" "${@}" >&2
}

function warning() {
	_color_log "WARNING" "${TTY_YELLOW}${@}${TTY_END}" >&2
}

function error() {
	_color_log "ERROR" "${TTY_RED}${@}${TTY_END}" >&2
}

function log_command() {
	info "  $ ${@}"
	"${@}"
}
