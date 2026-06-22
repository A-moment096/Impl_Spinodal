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
usage()   {
    echo -e "${CYAN}[USAGE]${RESET} $0 <source.c> <out-directory> [debug|release]"
    exit 1
}

# ── Usage ──────────────────────────────────────────────────────────────────

if [ $# -lt 2 ] || [ $# -gt 3 ]; then
    usage
fi

SRC="$1"
if [ ! -f "$SRC" ]; then
    error "$SRC is not a valid source file, please check the input source file." 
    usage 
fi

OUTDIR="$2"
if [ ! -d "$OUTDIR" ]; then
    error "$OUTDIR is not a valid directory, please check your output directory."
    usage
fi

BUILD_TYPE="${3:-debug}"
if [ -n "$build_type" ] && [ "$build_type" != "debug" ] && [ "$build_type" != "release" ]; then
    error "Third argument, if provided, must be 'debug' or 'release'"
    usage
fi

BASENAME=$(basename "$SRC" .c)

if [[ "$BASENAME" == *"v3"* ]]; then LIB="-lfftw3"; else LIB="-lmyfft"; fi
if [[ "$BUILD_TYPE" == "release" ]]; then FLAGS="-O2 -DNDEBUG"; else FLAGS="-g"; fi

mkdir -p "$OUTDIR"

C_OPTION=(-fdiagnostics-color=always -Wall -Wextra -pedantic $FLAGS \
    "$SRC" -o "$OUTDIR/$BASENAME" \
    -L"$(dirname "$SRC")/lib/linux" $LIB -lm)

info  "-------------------------"
info  "- Building information: -"
info  "-------------------------"
info  "Source:     $BASENAME.c"
info  "Output:     $OUTDIR"
info  "Build type: ${BOLD}$BUILD_TYPE${RESET}"
info  "Library:    $LIB"
info  "Build Opts: ${C_OPTION[@]}"
echo -e ""

gcc -fdiagnostics-color=always -Wall -Wextra -pedantic $FLAGS \
    "$SRC" -o "$OUTDIR/$BASENAME" \
    -L"$(dirname "$SRC")/lib/linux" $LIB -lm

if [[ $? -eq 0 ]]; then
    success "Build succeeded: $OUTDIR/$BASENAME"
else
    error "Build failed"
    exit 1
fi