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
SRC="$1"
OUTDIR="$2"
BUILD_TYPE="${3:-debug}"
BASENAME=$(basename "$SRC" .c)

if [[ "$BASENAME" == *"v3"* ]]; then LIB="-lfftw3"; else LIB="-lmyfft"; fi
if [[ "$BUILD_TYPE" == "release" ]]; then FLAGS="-O2 -DNDEBUG"; else FLAGS="-g"; fi

info  "-------------------------"
info  "- Building information: -"
info  "-------------------------"
info  "Source:     $BASENAME.c"
info  "Output:     $OUTDIR"
info  "Build type: ${BOLD}$BUILD_TYPE${RESET}"
info  "Library:    $LIB"
echo -e ""

mkdir -p "$OUTDIR"
gcc -fdiagnostics-color=always -Wall -Wextra -pedantic $FLAGS \
    "$SRC" -o "$OUTDIR/$BASENAME" \
    -L"$(dirname "$SRC")/lib/linux" $LIB -lm

if [[ $? -eq 0 ]]; then
    success "Build succeeded: $OUTDIR/$BASENAME"
else
    error "Build failed"
    exit 1
fi