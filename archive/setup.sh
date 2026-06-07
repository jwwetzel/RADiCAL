#!/usr/bin/env bash
# =============================================================================
# setup.sh — install prerequisites for the RADiCAL analysis workflow
# =============================================================================
#
# Installs:
#   1. Homebrew  (macOS package manager, if not already present)
#   2. ROOT      (CERN data analysis framework — required to run the macros)
#   3. Python scientific stack (optional — numpy, matplotlib, scipy, uproot)
#
# Usage (from anywhere):
#   bash Analysis/setup.sh            # interactive, asks before modifying shell
#   bash Analysis/setup.sh --no-python  # skip Python packages
#
# After installation, ROOT will be available in new terminal sessions.
# To activate it in the current session immediately:
#   source $(brew --prefix root)/bin/thisroot.sh
#
# Cluster usage (Iowa Argon / Lxplus):
#   See the "CLUSTER NOTES" section at the bottom of this script.
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
info()    { echo -e "${GREEN}[setup]${NC} $*"; }
warn()    { echo -e "${YELLOW}[setup]${NC} $*"; }
error()   { echo -e "${RED}[setup]${NC} $*" >&2; }
ask_yes() {                        # ask_yes "message" → returns 0 for yes
    read -r -p "$(echo -e "${YELLOW}[setup]${NC} $* [y/N] ")" ans
    [[ "$ans" =~ ^[Yy]$ ]]
}

INSTALL_PYTHON=true
for arg in "$@"; do
    [[ "$arg" == "--no-python" ]] && INSTALL_PYTHON=false
done

echo ""
echo "======================================================================"
echo "  RADiCAL analysis — prerequisite installer"
echo "  macOS $(sw_vers -productVersion)  $(uname -m)"
echo "======================================================================"
echo ""

# ---------------------------------------------------------------------------
# 1. Xcode Command Line Tools (required by Homebrew)
# ---------------------------------------------------------------------------
if ! xcode-select -p &>/dev/null; then
    info "Installing Xcode Command Line Tools (required by Homebrew)..."
    xcode-select --install
    echo ""
    warn "A dialog appeared — click 'Install' and wait for it to finish."
    warn "Then re-run this script."
    exit 0
fi
info "Xcode Command Line Tools: OK"

# ---------------------------------------------------------------------------
# 2. Homebrew
# ---------------------------------------------------------------------------
BREW=""
if command -v brew &>/dev/null; then
    BREW="$(command -v brew)"
    info "Homebrew already installed: $BREW"
elif [[ -x /opt/homebrew/bin/brew ]]; then
    BREW=/opt/homebrew/bin/brew
    info "Homebrew found at /opt/homebrew (not in PATH yet)"
elif [[ -x /usr/local/bin/brew ]]; then
    BREW=/usr/local/bin/brew
    info "Homebrew found at /usr/local (not in PATH yet)"
else
    info "Homebrew not found — installing now."
    echo ""
    warn "This will run the official Homebrew install script from:"
    warn "  https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh"
    echo ""
    if ask_yes "Proceed with Homebrew installation?"; then
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
        # Pick up the new brew binary
        if [[ -x /opt/homebrew/bin/brew ]]; then
            BREW=/opt/homebrew/bin/brew
        elif [[ -x /usr/local/bin/brew ]]; then
            BREW=/usr/local/bin/brew
        else
            error "Homebrew installation appears to have failed — brew binary not found."
            exit 1
        fi
        info "Homebrew installed at: $BREW"
    else
        error "Homebrew is required. Exiting."
        exit 1
    fi
fi

# Activate Homebrew in this shell session
eval "$($BREW shellenv)"

# ---------------------------------------------------------------------------
# 3. ROOT
# ---------------------------------------------------------------------------
info "Checking for ROOT..."
if $BREW list root &>/dev/null; then
    info "ROOT already installed via Homebrew."
elif command -v root &>/dev/null; then
    warn "ROOT found in PATH (not via Homebrew) — skipping installation."
    warn "  $(command -v root)"
else
    info "Installing ROOT via Homebrew (this downloads ~500 MB, may take 10–20 min)..."
    $BREW install root
    info "ROOT installed."
fi

ROOT_PREFIX=$($BREW --prefix root 2>/dev/null || true)
if [[ -z "$ROOT_PREFIX" ]]; then
    # ROOT might be keg-only or in an unusual location
    ROOT_PREFIX=$($BREW --cellar)/root/$($BREW list --versions root | awk '{print $2}')
fi
THISROOT="$ROOT_PREFIX/bin/thisroot.sh"

if [[ ! -f "$THISROOT" ]]; then
    # Fallback: find thisroot.sh
    THISROOT=$(find "$($BREW --cellar)/root" -name "thisroot.sh" 2>/dev/null | head -1 || true)
fi

if [[ -f "$THISROOT" ]]; then
    info "ROOT environment script: $THISROOT"
    # Activate ROOT in this session
    # shellcheck disable=SC1090
    source "$THISROOT"
    info "ROOT $(root --version 2>&1 | head -1) — active in this session."
else
    warn "Could not find thisroot.sh — you may need to source it manually."
fi

# ---------------------------------------------------------------------------
# 4. Python scientific packages (optional) — installed into a venv
# ---------------------------------------------------------------------------
#
# Homebrew Python 3.12+ is "externally managed" (PEP 668) and blocks
# system-wide pip installs.  We create a project-local virtual environment
# at .venv/ in the repository root instead.
# ---------------------------------------------------------------------------

# Resolve the repo root (one level up from Analysis/)
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VENV_DIR="$REPO_ROOT/.venv"

if $INSTALL_PYTHON; then
    PYTHON=""
    if command -v python3 &>/dev/null; then
        PYTHON=$(command -v python3)
    fi

    if [[ -z "$PYTHON" ]]; then
        warn "python3 not found — skipping Python packages."
        warn "Run 'brew install python' and re-run this script to add them."
    else
        info "Python: $PYTHON ($($PYTHON --version 2>&1))"

        # Create venv if it doesn't exist yet
        if [[ ! -d "$VENV_DIR" ]]; then
            info "Creating virtual environment at $VENV_DIR ..."
            "$PYTHON" -m venv "$VENV_DIR"
        else
            info "Virtual environment already exists at $VENV_DIR"
        fi

        VENV_PY="$VENV_DIR/bin/python3"
        VENV_PIP="$VENV_DIR/bin/pip"

        info "Installing into venv: numpy  matplotlib  scipy  uproot  awkward  hist"
        "$VENV_PIP" install --upgrade pip --quiet
        "$VENV_PIP" install --upgrade \
            numpy matplotlib scipy \
            uproot awkward hist \
            --quiet

        info "Python packages installed."
        "$VENV_PY" -c "import uproot;    print(f'  uproot     {uproot.__version__}')"
        "$VENV_PY" -c "import numpy;     print(f'  numpy      {numpy.__version__}')"
        "$VENV_PY" -c "import matplotlib; print(f'  matplotlib {matplotlib.__version__}')"
        "$VENV_PY" -c "import scipy;     print(f'  scipy      {scipy.__version__}')"
    fi
else
    info "Skipping Python packages (--no-python)."
fi

# ---------------------------------------------------------------------------
# 5. Shell configuration
# ---------------------------------------------------------------------------
echo ""
info "Configuring your shell so ROOT is available in all future sessions..."

SHELL_RC=""
if [[ "$SHELL" == */zsh ]] || [[ -f ~/.zshrc ]]; then
    SHELL_RC="$HOME/.zshrc"
elif [[ "$SHELL" == */bash ]]; then
    SHELL_RC="$HOME/.bash_profile"
fi

BREW_INIT_LINE="eval \"\$($BREW shellenv)\""
ROOT_INIT_LINE="[[ -f \"$THISROOT\" ]] && source \"$THISROOT\""

if [[ -n "$SHELL_RC" ]]; then
    # Add Homebrew to shell
    if ! grep -qF "$BREW_INIT_LINE" "$SHELL_RC" 2>/dev/null; then
        echo "" >> "$SHELL_RC"
        echo "# Homebrew" >> "$SHELL_RC"
        echo "$BREW_INIT_LINE" >> "$SHELL_RC"
        info "Added Homebrew init to $SHELL_RC"
    else
        info "Homebrew already configured in $SHELL_RC"
    fi
    # Add ROOT to shell
    if [[ -f "$THISROOT" ]] && ! grep -qF "thisroot.sh" "$SHELL_RC" 2>/dev/null; then
        echo "" >> "$SHELL_RC"
        echo "# ROOT (CERN data analysis framework)" >> "$SHELL_RC"
        echo "$ROOT_INIT_LINE" >> "$SHELL_RC"
        info "Added ROOT init to $SHELL_RC"
    elif grep -qF "thisroot.sh" "$SHELL_RC" 2>/dev/null; then
        info "ROOT already configured in $SHELL_RC"
    fi
else
    warn "Could not detect shell config file — add these lines manually:"
    echo ""
    echo "    $BREW_INIT_LINE"
    echo "    $ROOT_INIT_LINE"
    echo ""
fi

# Print venv activation reminder (not added to .zshrc automatically —
# the venv should be activated explicitly when doing Python work)
if [[ -d "$VENV_DIR" ]]; then
    echo ""
    info "Python venv created at: $VENV_DIR"
    info "Activate it when doing Python-based analysis:"
    echo ""
    echo "    source $VENV_DIR/bin/activate"
    echo ""
    info "To deactivate: type 'deactivate'"
fi

# ---------------------------------------------------------------------------
# 6. Quick sanity check
# ---------------------------------------------------------------------------
echo ""
echo "======================================================================"
info "Sanity check"
echo "======================================================================"
echo ""
if command -v root &>/dev/null; then
    root --version 2>&1 | head -1 | sed 's/^/  /'
    info "root binary: $(command -v root)"
    info "All good — you can now run the analysis:"
    echo ""
    echo "    bash Analysis/runAll.sh"
    echo ""
else
    warn "ROOT is installed but not yet in PATH for this shell session."
    warn "Either open a new terminal tab, or run:"
    echo ""
    if [[ -f "$THISROOT" ]]; then
        echo "    source \"$THISROOT\""
    fi
    echo ""
fi

# =============================================================================
# CLUSTER NOTES (not executed — read-only reference)
# =============================================================================
# Iowa Argon Cluster:
#   module load root/6.XX.XX       (check: module avail root)
#   OR: source /Shared/lss_yonel/.../thisroot.sh
#
# Lxplus (CERN):
#   LCG release:
#     source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh
#   Standalone ROOT:
#     source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.00/\
#            x86_64-el9-gcc13-opt/bin/thisroot.sh
#
# After sourcing, run from the repository root:
#   bash Analysis/runAll.sh
# =============================================================================
