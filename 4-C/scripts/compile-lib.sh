#!/bin/bash

# ── Colour definitions ─────────────────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
RESET='\033[0m'

# ── Helper functions ───────────────────────────────────────────────────────
info()    { echo -e "${CYAN}[INFO]${RESET}  $*"; }
success() { echo -e "${GREEN}[OK]${RESET}    $*"; }
warn()    { echo -e "${YELLOW}[WARN]${RESET}  $*"; }
error()   { echo -e "${RED}[ERROR]${RESET} $*" >&2; }

# ── Usage ──────────────────────────────────────────────────────────────────

SCRIPT_PATH="${BASH_SOURCE[0]}"
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"
WORKING_DIR=$(realpath $SCRIPT_DIR/..)
BUILD_DIR="$WORKING_DIR/build/linux"
LIB_DIR="$WORKING_DIR/lib/linux"

mkdir -p $BUILD_DIR
mkdir -p $LIB_DIR

C_OPTIONS=(-c -O2 -Wall -Wextra $WORKING_DIR/C_my_fft.c -o $BUILD_DIR/C_my_fft.o)

info "Compile options: ${C_OPTIONS[@]}"

gcc "${C_OPTIONS[@]}"
if [[ $? -ne 0 ]]; then
    error "Compile failed with exit code $?"
    exit 1
fi

ar rcs $LIB_DIR/libmyfft.a $BUILD_DIR/C_my_fft.o
if [[ $? -eq 0 ]]; then
    success "Compile succeeded: $LIB_DIR/libmyfft.a"
else
    error "Compile failed with exit code $?"
    exit 1
fi