#!/usr/bin/env bash
# Ribologic setup
#
# Supported:
#   - macOS (Intel and Apple Silicon)
#   - Linux (x86_64 and arm64)
#   - Native Windows through MSYS2 MinGW64
#   - WSL2 (treated as Linux)
#
# Native Windows prerequisite:
#   1. Install MSYS2: https://www.msys2.org/
#   2. Open "MSYS2 MinGW x64"
#   3. Clone this repository and run:
#        ./setup.sh

set -euo pipefail

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

ok()   { echo -e "${GREEN}✓ $*${NC}"; }
warn() { echo -e "${YELLOW}! $*${NC}"; }
die()  { echo -e "${RED}✗ $*${NC}" >&2; exit 1; }

# ViennaRNA version for macOS/Linux source builds.
VIENNARNA_VERSION="${VIENNARNA_VERSION:-2.6.4}"

# ViennaRNA version for Windows MinGW source builds.
VIENNARNA_WINDOWS_VERSION="${VIENNARNA_WINDOWS_VERSION:-2.7.0}"

# Always work from the repository root, where this script lives.
cd "$(dirname "${BASH_SOURCE[0]}")"

PROJECT_DIR="$(pwd)"
VENDOR_DIR="$PROJECT_DIR/vendor/RNAlib"

echo "Ribologic Setup"
echo "================"
echo ""

# ── Platform detection ──────────────────────────────────────────

OS_RAW="$(uname -s)"
ARCH_RAW="$(uname -m)"

case "$OS_RAW" in
    Darwin)
        OS="macos"
        ;;
    Linux)
        OS="linux"
        ;;
    MINGW*|MSYS*)
        OS="windows"
        ;;
    CYGWIN*)
        die "Cygwin is not supported. Install MSYS2 from https://www.msys2.org/ and run this script from the 'MSYS2 MinGW x64' terminal."
        ;;
    *)
        die "Unsupported operating system: $OS_RAW"
        ;;
esac

case "$ARCH_RAW" in
    x86_64|amd64)
        ARCH="x86_64"
        ;;
    arm64|aarch64)
        ARCH="arm64"
        ;;
    *)
        die "Unsupported CPU architecture: $ARCH_RAW"
        ;;
esac

if [[ "$OS" == "windows" && "$ARCH" != "x86_64" ]]; then
    die "Native Windows setup currently supports x86_64 only."
fi

# Map the local machine to the exact Rust target/build.rs directory name.
case "$OS:$ARCH" in
    macos:x86_64)
        RUST_TARGET="x86_64-apple-darwin"
        ;;
    macos:arm64)
        RUST_TARGET="aarch64-apple-darwin"
        ;;
    linux:x86_64)
        RUST_TARGET="x86_64-unknown-linux-gnu"
        ;;
    linux:arm64)
        RUST_TARGET="aarch64-unknown-linux-gnu"
        ;;
    windows:x86_64)
        RUST_TARGET="x86_64-pc-windows-gnu"
        ;;
    *)
        die "No Rust target mapping is available for $OS ($ARCH)."
        ;;
esac

# This exact location must match what build.rs expects.
PREBUILT_DIR="$VENDOR_DIR/prebuilt/$RUST_TARGET"

# Names used by setup functions below.
NATIVE_VENDOR_LIB_DIR="$PREBUILT_DIR"
WINDOWS_TARGET="$RUST_TARGET"
WINDOWS_VENDOR_LIB_DIR="$PREBUILT_DIR"

ok "Detected $OS ($ARCH)"
ok "Rust target: $RUST_TARGET"
ok "Native library directory: $PREBUILT_DIR"

# ── Common helpers ──────────────────────────────────────────────

have() {
    command -v "$1" >/dev/null 2>&1
}

need_curl() {
    have curl || die "curl is required but was not found."
}

cpu_count() {
    if have nproc; then
        nproc
    elif have sysctl; then
        sysctl -n hw.ncpu
    else
        echo 2
    fi
}

SUDO=""
if [[ "$OS" != "windows" && "$(id -u)" -ne 0 ]]; then
    if have sudo; then
        SUDO="sudo"
    else
        warn "sudo was not found. System package installation may fail."
    fi
fi

# Make the Windows GitHub Actions Rust installation visible from MSYS2.
#
# On GitHub-hosted Windows runners, USERPROFILE is generally:
# C:\Users\runneradmin
#
# cygpath converts it into the MSYS2 form:
# /c/Users/runneradmin
if [[ "$OS" == "windows" ]]; then
    have cygpath ||
        die "cygpath was not found. Run this script from MSYS2 MINGW64."

    WINDOWS_HOME="$(cygpath -u "$USERPROFILE")"

    export CARGO_HOME="${CARGO_HOME:-$WINDOWS_HOME/.cargo}"
    export RUSTUP_HOME="${RUSTUP_HOME:-$WINDOWS_HOME/.rustup}"
    export PATH="$CARGO_HOME/bin:$PATH"

    export LIBCLANG_PATH="${LIBCLANG_PATH:-/mingw64/bin}"

    # Persist these variables for subsequent GitHub Actions steps.
    #
    # GITHUB_PATH is processed by GitHub Actions after setup.sh's step ends.
    # cygpath -w makes the path readable by the Windows runner.
    if [[ -n "${GITHUB_PATH:-}" ]]; then
        cygpath -w "$CARGO_HOME/bin" >> "$GITHUB_PATH"
        cygpath -w "/mingw64/bin" >> "$GITHUB_PATH"
    fi

    if [[ -n "${GITHUB_ENV:-}" ]]; then
        echo "LIBCLANG_PATH=/mingw64/bin" >> "$GITHUB_ENV"
        echo "CARGO_HOME=$(cygpath -w "$CARGO_HOME")" >> "$GITHUB_ENV"
        echo "RUSTUP_HOME=$(cygpath -w "$RUSTUP_HOME")" >> "$GITHUB_ENV"
    fi
fi





# ── Rust setup ──────────────────────────────────────────────────

setup_rust() {
    if [[ "$OS" == "windows" ]]; then
        have cygpath ||
            die "cygpath was not found. Run this script from MSYS2 MinGW x64."

        # If Cargo is already visible, preserve the existing configuration.
        #
        # This supports:
        # - Rust installed through rustup on Windows
        # - a GitHub Actions Windows runner
        # - Cargo already configured by the developer
        if have cargo; then
            CARGO_BIN_DIR="$(dirname "$(command -v cargo)")"
            export PATH="$CARGO_BIN_DIR:$PATH"

            # Only infer CARGO_HOME if the user did not explicitly set it.
            if [[ -z "${CARGO_HOME:-}" ]]; then
                export CARGO_HOME="$(dirname "$CARGO_BIN_DIR")"
            fi
        else
            # Standard rustup location for the current Windows user.
            #
            # USERPROFILE example:
            # C:\Users\Alice
            #
            # Converted MSYS2 path:
            # /c/Users/Alice
            WINDOWS_HOME="$(cygpath -u "${USERPROFILE:-}")"

            [[ -n "$WINDOWS_HOME" && "$WINDOWS_HOME" != "." ]] ||
                die "USERPROFILE is not set; cannot determine the Windows Rust installation directory."

            export CARGO_HOME="${CARGO_HOME:-$WINDOWS_HOME/.cargo}"
            export RUSTUP_HOME="${RUSTUP_HOME:-$WINDOWS_HOME/.rustup}"
            export PATH="$CARGO_HOME/bin:$PATH"
        fi
    else
        export CARGO_HOME="${CARGO_HOME:-$HOME/.cargo}"
        export RUSTUP_HOME="${RUSTUP_HOME:-$HOME/.rustup}"
        export PATH="$CARGO_HOME/bin:$PATH"
    fi

    if ! have cargo && [[ -f "$CARGO_HOME/env" ]]; then
        # shellcheck disable=SC1090
        source "$CARGO_HOME/env"
    fi

    if ! have cargo; then
        need_curl

        echo "Cargo was not found; installing Rust..."

        curl --proto '=https' --tlsv1.2 -sSf \
            https://sh.rustup.rs | sh -s -- -y

        export PATH="$CARGO_HOME/bin:$PATH"

        if [[ -f "$CARGO_HOME/env" ]]; then
            # shellcheck disable=SC1090
            source "$CARGO_HOME/env"
        fi
    fi

    have cargo ||
        die "Rust installation failed: cargo was not found."

    have rustup ||
        die "Rust installation failed: rustup was not found."

    if [[ "$OS" == "windows" ]]; then
        echo "Installing/selecting Windows GNU Rust toolchain..."

        rustup toolchain install stable-x86_64-pc-windows-gnu
        rustup override set stable-x86_64-pc-windows-gnu

        rustup target add \
            --toolchain stable-x86_64-pc-windows-gnu \
            x86_64-pc-windows-gnu

        # Only relevant in GitHub Actions; harmless for local users.
        if [[ -n "${GITHUB_PATH:-}" ]]; then
            cygpath -w "$(dirname "$(command -v cargo)")" >> "$GITHUB_PATH"
            cygpath -w "/mingw64/bin" >> "$GITHUB_PATH"
        fi
    else
        rustup target add "$RUST_TARGET"
    fi

    echo "Cargo location:  $(command -v cargo)"
    echo "Rustup location: $(command -v rustup)"

    ok "Rust ready: $(cargo --version)"
    ok "Rust compiler: $(rustc --version)"
}





# ── Library checks ──────────────────────────────────────────────

REQUIRED_LIBRARIES=(
    libRNA.a
    libgsl.a
    libgslcblas.a
    libmpfr.a
    libgmp.a
)

native_vendor_complete() {
    local lib

    for lib in "${REQUIRED_LIBRARIES[@]}"; do
        [[ -f "$NATIVE_VENDOR_LIB_DIR/$lib" ]] || return 1
    done

    return 0
}

windows_vendor_complete() {
    local lib

    for lib in "${REQUIRED_LIBRARIES[@]}"; do
        [[ -f "$WINDOWS_VENDOR_LIB_DIR/$lib" ]] || return 1
    done

    return 0
}

# ── macOS ───────────────────────────────────────────────────────

copy_brew_lib() {
    local formula="$1"
    local lib="$2"
    local prefix
    local source_file

    prefix="$(brew --prefix "$formula")"
    source_file="$prefix/lib/$lib"

    # Homebrew occasionally places archives elsewhere in the formula cellar.
    if [[ ! -f "$source_file" ]]; then
        source_file="$(
            find "$(brew --cellar "$formula")" \
                -type f \
                -name "$lib" \
                -print \
                -quit 2>/dev/null || true
        )"
    fi

    if [[ -n "$source_file" && -f "$source_file" ]]; then
        cp "$source_file" "$NATIVE_VENDOR_LIB_DIR/"
        return 0
    fi

    return 1
}

setup_macos() {
    if ! xcode-select -p >/dev/null 2>&1; then
        echo "Installing Xcode Command Line Tools..."
        xcode-select --install || true
        die "Finish installing Xcode Command Line Tools, then re-run ./setup.sh."
    fi

    # Find an existing Homebrew installation if brew is not already on PATH.
    if ! have brew; then
        for candidate in /opt/homebrew/bin/brew /usr/local/bin/brew; do
            if [[ -x "$candidate" ]]; then
                eval "$("$candidate" shellenv)"
                break
            fi
        done
    fi

    # Install Homebrew if needed.
    if ! have brew; then
        need_curl

        echo "Installing Homebrew..."
        /bin/bash -c \
            "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

        if [[ -x /opt/homebrew/bin/brew ]]; then
            eval "$(/opt/homebrew/bin/brew shellenv)"
        elif [[ -x /usr/local/bin/brew ]]; then
            eval "$(/usr/local/bin/brew shellenv)"
        fi
    fi

    have brew || die "Homebrew installation failed."

    echo "Installing macOS dependencies..."

    # Do NOT add `viennarna` here: Homebrew no longer has that formula.
    brew install \
        gsl \
        mpfr \
        gmp \
        llvm \
        autoconf \
        automake \
        libtool \
        pkg-config \
        bison \
        flex

    mkdir -p "$NATIVE_VENDOR_LIB_DIR" "$VENDOR_DIR/include"

    # Make pkg-config and the compiler find Homebrew packages.
    export PATH="$(brew --prefix llvm)/bin:$PATH"
    export PKG_CONFIG_PATH="$(brew --prefix gsl)/lib/pkgconfig:$(brew --prefix mpfr)/lib/pkgconfig:$(brew --prefix gmp)/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
    export CPPFLAGS="-I$(brew --prefix gsl)/include -I$(brew --prefix mpfr)/include -I$(brew --prefix gmp)/include ${CPPFLAGS:-}"
    export LDFLAGS="-L$(brew --prefix gsl)/lib -L$(brew --prefix mpfr)/lib -L$(brew --prefix gmp)/lib ${LDFLAGS:-}"

    echo "Copying GSL, MPFR, and GMP static libraries..."

    for pair in \
        "gsl:libgsl.a" \
        "gsl:libgslcblas.a" \
        "mpfr:libmpfr.a" \
        "gmp:libgmp.a"
    do
        local formula="${pair%%:*}"
        local library="${pair##*:}"

        copy_brew_lib "$formula" "$library" ||
            die "Could not find static library $library from Homebrew formula $formula."
    done

    # ViennaRNA must be built from source on macOS.
    build_viennarna_native

    ok "macOS native libraries ready"
}


# ── Linux ───────────────────────────────────────────────────────

install_linux_packages() {
    if have apt-get; then
        echo "Using apt..."

        $SUDO apt-get update
        $SUDO apt-get install -y \
            build-essential \
            curl \
            pkg-config \
            clang \
            libclang-dev \
            autoconf \
            automake \
            libtool \
            bison \
            flex \
            libgsl-dev \
            libmpfr-dev \
            libgmp-dev

    elif have dnf; then
        echo "Using dnf..."

        $SUDO dnf install -y \
            gcc \
            gcc-c++ \
            make \
            curl \
            pkgconf-pkg-config \
            clang \
            clang-libs \
            autoconf \
            automake \
            libtool \
            bison \
            flex \
            gsl-devel \
            gsl-static \
            mpfr-devel \
            gmp-devel \
            gmp-static

    elif have yum; then
        echo "Using yum..."

        $SUDO yum install -y \
            gcc \
            gcc-c++ \
            make \
            curl \
            pkgconfig \
            clang \
            autoconf \
            automake \
            libtool \
            bison \
            flex \
            gsl-devel \
            mpfr-devel \
            gmp-devel

    elif have pacman; then
        echo "Using pacman..."

        $SUDO pacman -Sy --needed --noconfirm \
            base-devel \
            curl \
            pkgconf \
            clang \
            autoconf \
            automake \
            libtool \
            bison \
            flex \
            gsl \
            mpfr \
            gmp

    elif have zypper; then
        echo "Using zypper..."

        $SUDO zypper --non-interactive install \
            gcc \
            gcc-c++ \
            make \
            curl \
            pkg-config \
            clang \
            libclang-devel \
            autoconf \
            automake \
            libtool \
            bison \
            flex \
            gsl-devel \
            mpfr-devel \
            gmp-devel

    else
        die "No supported Linux package manager found. Install a compiler, make, curl, libclang, GSL, MPFR, GMP, and autotools manually."
    fi
}

find_static_lib() {
    local name="$1"
    local path

    path="$(gcc -print-file-name="$name" 2>/dev/null || true)"

    if [[ -n "$path" && -f "$path" ]]; then
        echo "$path"
        return 0
    fi

    for dir in /usr/lib64 /usr/lib /usr/local/lib /usr/lib/*-linux-gnu; do
        if [[ -f "$dir/$name" ]]; then
            echo "$dir/$name"
            return 0
        fi
    done

    return 1
}

setup_linux() {
    install_linux_packages

    mkdir -p "$NATIVE_VENDOR_LIB_DIR" "$VENDOR_DIR/include"

    # Most Linux distributions do not package static ViennaRNA.
    if [[ ! -f "$NATIVE_VENDOR_LIB_DIR/libRNA.a" ]]; then
        build_viennarna_native
    fi

    echo "Copying Linux GSL, MPFR, and GMP static libraries..."

    local lib
    local source_file

    for lib in libgsl.a libgslcblas.a libmpfr.a libgmp.a; do
        source_file="$(find_static_lib "$lib")" ||
            die "Could not find $lib. Install the relevant static development package."

        cp "$source_file" "$NATIVE_VENDOR_LIB_DIR/"
    done

    ok "Linux native libraries ready"
}

# ── Build ViennaRNA for macOS/Linux ─────────────────────────────

build_viennarna_native() {
    need_curl

    local temp_dir
    local rna_archive

    temp_dir="$(mktemp -d)"

    echo "Building ViennaRNA $VIENNARNA_VERSION from source..."

    curl --fail --location \
        --retry 5 \
        --retry-delay 5 \
        -o "$temp_dir/viennarna.tar.gz" \
        "https://github.com/ViennaRNA/ViennaRNA/releases/download/v${VIENNARNA_VERSION}/ViennaRNA-${VIENNARNA_VERSION}.tar.gz"

    tar -xzf "$temp_dir/viennarna.tar.gz" -C "$temp_dir"

    pushd "$temp_dir/ViennaRNA-${VIENNARNA_VERSION}" >/dev/null

    if [[ "$OS" == "macos" ]]; then
        export CFLAGS="${CFLAGS:-} -Wno-error -Wno-implicit-function-declaration -Wno-int-conversion"
    fi

    ./configure \
        --prefix="$NATIVE_VENDOR_LIB_DIR" \
        --disable-shared \
        --enable-static \
        --with-pic \
        --without-perl \
        --without-python \
        --without-doc \
        --without-kinfold \
        --without-forester \
        --without-rnalocmin \
        --disable-lto \
        --disable-openmp

    make -j"$(cpu_count)"

    # Do not use `make install`.
    #
    # ViennaRNA's documentation installation target can fail on macOS.
    # We only need the static archive and C headers for the Rust build.
    rna_archive="$(find . -type f -name libRNA.a -print -quit)"

    [[ -n "$rna_archive" && -f "$rna_archive" ]] ||
        die "ViennaRNA build completed but did not produce libRNA.a."

    mkdir -p "$NATIVE_VENDOR_LIB_DIR" "$VENDOR_DIR/include"

    # build.rs expects the archive directly in the target-specific directory.
    cp "$rna_archive" "$NATIVE_VENDOR_LIB_DIR/libRNA.a"

    # Preserve the expected include layout:
    # vendor/RNAlib/include/ViennaRNA/*.h
    #
    # The configured source tree includes generated headers needed by ViennaRNA.
    rm -rf "$VENDOR_DIR/include/ViennaRNA"
    cp -R "$PWD/src/ViennaRNA" "$VENDOR_DIR/include/"

    popd >/dev/null
    rm -rf "$temp_dir"

    [[ -f "$NATIVE_VENDOR_LIB_DIR/libRNA.a" ]] ||
        die "ViennaRNA library was not copied to $NATIVE_VENDOR_LIB_DIR/libRNA.a"

    [[ -d "$VENDOR_DIR/include/ViennaRNA" ]] ||
        die "ViennaRNA headers were not copied to $VENDOR_DIR/include/ViennaRNA"

    ok "ViennaRNA static library and headers staged"
}


# ── Native Windows through MSYS2 MinGW64 ────────────────────────

setup_windows() {
    [[ "${MSYSTEM:-}" == "MINGW64" ]] ||
        die "Run this script from the 'MSYS2 MinGW x64' terminal, not PowerShell, CMD, Git Bash, WSL, or the plain MSYS2 terminal."

    export PATH="/mingw64/bin:$PATH"
    export PKG_CONFIG_PATH="/mingw64/lib/pkgconfig"

    export LIBCLANG_PATH="/mingw64/bin"

    if [[ -n "${GITHUB_ENV:-}" ]]; then
        echo "LIBCLANG_PATH=/mingw64/bin" >> "$GITHUB_ENV"
    fi


    echo "Installing MSYS2 MinGW dependencies..."

    pacman -Sy --needed --noconfirm \
        mingw-w64-x86_64-toolchain \
        mingw-w64-x86_64-clang \
        mingw-w64-x86_64-pkgconf \
        mingw-w64-x86_64-gmp \
        mingw-w64-x86_64-mpfr \
        mingw-w64-x86_64-gsl \
        git \
        curl \
        make \
        autoconf \
        automake \
        libtool \
        bison \
        flex

    #rustup target add "$WINDOWS_TARGET"

    local temp_dir
    local rna_archive

    temp_dir="$(mktemp -d)"

    echo "Building ViennaRNA $VIENNARNA_WINDOWS_VERSION for Windows MinGW..."

    pushd "$temp_dir" >/dev/null

    curl --fail --location \
        --retry 5 \
        --retry-delay 5 \
        --connect-timeout 30 \
        --max-time 300 \
        -o "ViennaRNA-${VIENNARNA_WINDOWS_VERSION}.tar.gz" \
        "https://github.com/ViennaRNA/ViennaRNA/releases/download/v${VIENNARNA_WINDOWS_VERSION}/ViennaRNA-${VIENNARNA_WINDOWS_VERSION}.tar.gz"

    tar -xf "ViennaRNA-${VIENNARNA_WINDOWS_VERSION}.tar.gz"
    cd "ViennaRNA-${VIENNARNA_WINDOWS_VERSION}"

    export CC=gcc
    export CXX=g++
    export CPPFLAGS="-I$PWD/src -I/mingw64/include"
    export CFLAGS="-std=gnu17 -fno-lto"
    export CXXFLAGS="-std=gnu++17 -fno-lto -include cstdint"
    export LDFLAGS="-fno-lto -L/mingw64/lib -lws2_32 -lwinmm"

    ./configure \
        --prefix="/mingw64" \
        --disable-shared \
        --enable-static \
        --without-openmp \
        --without-swig \
        --without-svm \
        --without-forester

    make -j1

    rna_archive="$(find . -type f -name libRNA.a -print -quit)"

    [[ -n "$rna_archive" ]] ||
        die "ViennaRNA build did not produce libRNA.a."

    mkdir -p "$WINDOWS_VENDOR_LIB_DIR"

    cp "$rna_archive" "$WINDOWS_VENDOR_LIB_DIR/libRNA.a"
    cp /mingw64/lib/libgmp.a "$WINDOWS_VENDOR_LIB_DIR/libgmp.a"
    cp /mingw64/lib/libmpfr.a "$WINDOWS_VENDOR_LIB_DIR/libmpfr.a"
    cp /mingw64/lib/libgsl.a "$WINDOWS_VENDOR_LIB_DIR/libgsl.a"
    cp /mingw64/lib/libgslcblas.a "$WINDOWS_VENDOR_LIB_DIR/libgslcblas.a"

    popd >/dev/null
    rm -rf "$temp_dir"

    windows_vendor_complete ||
        die "Windows libraries were not staged correctly in $WINDOWS_VENDOR_LIB_DIR."

    ok "Windows GNU libraries ready in $WINDOWS_VENDOR_LIB_DIR"
}

# ── Main ────────────────────────────────────────────────────────

setup_rust

case "$OS" in
    macos)
        if native_vendor_complete; then
            ok "macOS vendored libraries already present"
        else
            setup_macos
        fi

        native_vendor_complete ||
            die "macOS libraries were not staged correctly in $NATIVE_VENDOR_LIB_DIR."
        ;;

    linux)
        if native_vendor_complete; then
            ok "Linux vendored libraries already present"
        else
            setup_linux
        fi

        native_vendor_complete ||
            die "Linux libraries were not staged correctly in $NATIVE_VENDOR_LIB_DIR."
        ;;

    windows)
        if windows_vendor_complete; then
            ok "Windows GNU vendored libraries already present"
        else
            setup_windows
        fi

        windows_vendor_complete ||
            die "Windows libraries were not staged correctly in $WINDOWS_VENDOR_LIB_DIR."
        ;;
esac

echo ""
ok "Setup complete for $OS ($ARCH)"
echo ""

if [[ "$OS" == "windows" ]]; then
    echo "Run the program with:"
    echo "  cargo run --target $RUST_TARGET"
    echo ""
    echo "Run tests with:"
    echo "  cargo test --target $RUST_TARGET"
else
    echo "Run the program with:"
    echo "  cargo run"
    echo ""
    echo "Run tests with:"
    echo "  cargo test"
fi
