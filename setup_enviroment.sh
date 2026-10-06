#!/bin/bash

set -e

echo "Setting up Ribologic environment..."

# --------------------------------------------------
# Check Python version
# --------------------------------------------------

if ! command -v python3.12 &> /dev/null; then
    echo "ERROR: Python 3.12 is required."
    exit 1
fi

PYTHON=python3.12

echo "Python:"
$PYTHON --version


# --------------------------------------------------
# Detect operating system
# --------------------------------------------------

OS="$(uname -s)"
ARCH="$(uname -m)"

echo "Operating system: $OS"
echo "Architecture:     $ARCH"


# --------------------------------------------------
# Select requirements file
# --------------------------------------------------

case "$OS-$ARCH" in

    Darwin-arm64)
        REQUIREMENTS="requirements/macos-arm64.txt"
        ;;

    Darwin-x86_64)
        REQUIREMENTS="requirements/macos-x86_64.txt"
        ;;

    Linux-x86_64)
        REQUIREMENTS="requirements/linux-x86_64.txt"
        ;;

    Linux-aarch64)
        REQUIREMENTS="requirements/linux-aarch64.txt"
        ;;

    *)
        echo "ERROR: Unsupported platform: $OS-$ARCH"
        exit 1
        ;;

esac

echo "Using requirements:"
echo "$REQUIREMENTS"


# --------------------------------------------------
# Create virtual environment
# --------------------------------------------------

if [ ! -d ".venv" ]; then

    echo "Creating Python 3.12 virtual environment..."

    $PYTHON -m venv .venv

else

    echo ".venv already exists."

fi


# --------------------------------------------------
# Activate environment
# --------------------------------------------------

source .venv/bin/activate


# --------------------------------------------------
# Install dependencies
# --------------------------------------------------

echo "Installing dependencies..."

python -m pip install --upgrade pip

python -m pip install -r "$REQUIREMENTS"


# --------------------------------------------------
# Finished
# --------------------------------------------------

echo ""
echo "===================================="
echo "Ribologic environment ready!"
echo "===================================="

echo ""
echo "Python:"
python --version

echo ""
echo "Platform:"
echo "$OS-$ARCH"

echo ""
echo "Requirements:"
echo "$REQUIREMENTS"

echo ""
echo "Activate with:"
echo "source .venv/bin/activate"
