#!/usr/bin/env bash
# =============================================================================
# install.sh -- set up everything Shiny-Sturgeon needs on a Linux machine.
#
#   ./install.sh                 # interactive: asks before each big step
#   ./install.sh --yes --accept-sturgeon-license   # unattended
#
# What it does (each step can be skipped, see --help):
#   1. installs the system libraries R needs to build its packages (apt, sudo)
#   2. creates two conda environments from envs/:
#        sturgeon     samtools, modkit and Sturgeon v1 (Python 3.8)
#        sturgeon-v2  sturgeon-v2 (Python 3.11)
#      Sturgeon v1 needs Python older than 3.10 and sturgeon-v2 needs 3.10 or
#      newer, so they cannot share one environment.
#   3. installs the R packages at the versions in renv.lock (renv::restore())
#   4. downloads the Sturgeon models and the Sturgeon v1 reference files
#      (R_scripts/). These fall under Sturgeon's own license, so they are
#      downloaded from their original sources, not bundled with this app.
#   5. writes config_local.yml, so the app opens with every tool and reference
#      file filled in. Setup then only needs the MinKNOW folder, the barcode
#      and an output folder.
#
# Running it again is safe: existing environments are updated, files that
# are already downloaded are kept, and the previous config_local.yml is kept
# as config_local.yml.bak.
# =============================================================================
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$REPO_DIR/resources"
ASSUME_YES=0
ACCEPT_LICENSE=0
DO_SYSTEM=1
DO_CONDA=1
DO_R=1
DO_DOWNLOADS=1
DO_V2=1

# Pinned sources (the versions in the README's "Tested with" table)
STURGEON_V1_COMMIT="80c33fbdf304ee2516eb2602433c02150c254a78"
STURGEON_V1_RAW="https://raw.githubusercontent.com/UMCUGenetics/sturgeon/$STURGEON_V1_COMMIT"
STURGEON_V1_FILES=(probelocs_chm13.bed chm13_v2.0_bins.rds bins_export.bed color_translation.rds
                   relevant_genes_with_chm13v2_1Mb_bin_nrs.bed)
MODEL_V1_URL="https://www.dropbox.com/s/yzca4exl40x9ukw/general.zip?dl=1"
MODEL_V2_ZENODO_RECORD="22899108"
MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh"
# Zenodo refuses anonymous script requests from some networks ("unusual
# traffic", HTTP 403), but accepts them once the client names itself
USER_AGENT="Shiny-Sturgeon-installer"
ENV_V1="sturgeon"
ENV_V2="sturgeon-v2"

usage() {
  cat <<EOF
Usage: ./install.sh [options]

Options:
  --data-dir DIR             where models and reference files go (default: $REPO_DIR/resources)
  --yes, -y                  don't ask before each step
  --accept-sturgeon-license  accept the Sturgeon license without being asked (needed with --yes)
  --no-v2                    skip sturgeon-v2 (its environment and model)
  --skip-system              don't install system libraries (no sudo needed)
  --skip-conda               don't create or update the conda environments
  --skip-r                   don't install R packages
  --skip-downloads           don't download models and reference files
  -h, --help                 show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-dir) DATA_DIR="$2"; shift 2 ;;
    --data-dir=*) DATA_DIR="${1#*=}"; shift ;;
    -y|--yes) ASSUME_YES=1; shift ;;
    --accept-sturgeon-license) ACCEPT_LICENSE=1; shift ;;
    --no-v2) DO_V2=0; shift ;;
    --skip-system) DO_SYSTEM=0; shift ;;
    --skip-conda) DO_CONDA=0; shift ;;
    --skip-r) DO_R=0; shift ;;
    --skip-downloads) DO_DOWNLOADS=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

mkdir -p "$DATA_DIR"
DATA_DIR="$(cd "$DATA_DIR" && pwd)"

# ---- output helpers ---------------------------------------------------------
if [[ -t 1 ]]; then BOLD=$'\e[1m'; GREEN=$'\e[32m'; YELLOW=$'\e[33m'; RED=$'\e[31m'; RESET=$'\e[0m'
else BOLD=""; GREEN=""; YELLOW=""; RED=""; RESET=""; fi
step()  { echo; echo "${BOLD}==> $*${RESET}"; }
ok()    { echo "    ${GREEN}OK${RESET}  $*"; }
warn()  { echo "    ${YELLOW}WARNING${RESET}  $*"; WARNINGS+=("$*"); }
fail()  { echo "${RED}ERROR:${RESET} $*" >&2; exit 1; }
WARNINGS=()

ask() {  # ask "question" -> 0 for yes
  [[ $ASSUME_YES -eq 1 ]] && return 0
  local reply
  read -r -p "    $1 [Y/n] " reply || reply=""
  [[ -z "$reply" || "$reply" =~ ^[Yy] ]]
}

download() {  # download URL DEST -- skips files that are already there
  local url="$1" dest="$2"
  if [[ -s "$dest" ]]; then ok "already present: $dest"; return 0; fi
  echo "    downloading $(basename "$dest") ..."
  if curl -fL --retry 3 -A "$USER_AGENT" --progress-bar -o "$dest.part" "$url"; then
    mv "$dest.part" "$dest"
    ok "$dest"
  else
    rm -f "$dest.part"
    warn "could not download $url -> $dest"
    return 1
  fi
}

is_zip() { [[ -s "$1" ]] && head -c 2 "$1" | grep -q "PK"; }

[[ "$(uname -s)" == "Linux" ]] || fail "Shiny-Sturgeon and its tools (modkit, Sturgeon) run on Linux only."
command -v curl >/dev/null || fail "curl is needed. Install it with: sudo apt-get install curl"

echo "${BOLD}Shiny-Sturgeon installer${RESET}"
echo "  app folder:          $REPO_DIR"
echo "  models & references: $DATA_DIR"

# ---- 1. system libraries ----------------------------------------------------
SYSTEM_PACKAGES=(libcurl4-openssl-dev libssl-dev libxml2-dev libbz2-dev liblzma-dev zlib1g-dev libdeflate-dev)
if [[ $DO_SYSTEM -eq 1 ]]; then
  step "1/5  System libraries for building R packages"
  if command -v apt-get >/dev/null; then
    missing=()
    for p in "${SYSTEM_PACKAGES[@]}"; do
      dpkg -s "$p" >/dev/null 2>&1 || missing+=("$p")
    done
    if [[ ${#missing[@]} -eq 0 ]]; then
      ok "all installed"
    elif ask "Install ${missing[*]} with apt-get (needs your sudo password)?"; then
      SUDO=""; [[ $EUID -ne 0 ]] && SUDO="sudo"
      $SUDO apt-get update -qq
      $SUDO apt-get install -y "${missing[@]}"
      ok "installed"
    else
      warn "skipped system libraries; R packages may fail to build"
    fi
  else
    warn "no apt-get (not Ubuntu/Debian). Install the equivalents of: ${SYSTEM_PACKAGES[*]}"
  fi
fi

# ---- 2. conda environments --------------------------------------------------
find_conda() {
  command -v conda >/dev/null && { command -v conda; return; }
  for base in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3" /opt/conda; do
    [[ -x "$base/bin/conda" ]] && { echo "$base/bin/conda"; return; }
  done
}

CONDA="$(find_conda || true)"
if [[ $DO_CONDA -eq 1 ]]; then
  step "2/5  Conda environments ($ENV_V1$( [[ $DO_V2 -eq 1 ]] && echo ", $ENV_V2"))"
  if [[ -z "$CONDA" ]]; then
    echo "    conda was not found."
    if ask "Install Miniforge (a small conda) into $HOME/miniforge3?"; then
      installer="$(mktemp --suffix=.sh)"
      curl -fsSL -A "$USER_AGENT" -o "$installer" "$MINIFORGE_URL" || fail "could not download Miniforge from $MINIFORGE_URL"
      bash "$installer" -b -p "$HOME/miniforge3"
      rm -f "$installer"
      CONDA="$HOME/miniforge3/bin/conda"
      ok "Miniforge installed in $HOME/miniforge3"
    else
      fail "conda is needed for samtools, modkit and Sturgeon. Install Miniforge or rerun with --skip-conda."
    fi
  fi
  ok "using $CONDA"

  env_exists() { "$CONDA" env list | awk '{print $1}' | grep -qx "$1"; }
  make_env() {  # make_env NAME FILE
    if env_exists "$1"; then
      echo "    updating existing environment '$1' (this can take a while) ..."
      "$CONDA" env update -n "$1" -f "$2" --prune
    else
      echo "    creating environment '$1' (this can take a while) ..."
      "$CONDA" env create -n "$1" -f "$2"
    fi
    ok "environment '$1'"
  }
  make_env "$ENV_V1" "$REPO_DIR/envs/sturgeon.yml"
  [[ $DO_V2 -eq 1 ]] && make_env "$ENV_V2" "$REPO_DIR/envs/sturgeon-v2.yml"
fi

# Locate the environments' bin/ folders (also when --skip-conda was used)
env_bin() {
  [[ -n "$CONDA" ]] || return 0
  "$CONDA" env list | awk -v n="$1" '$1 == n {print $NF}' | head -1 | sed 's|$|/bin|'
}
BIN_V1="$(env_bin "$ENV_V1")"
BIN_V2=""; [[ $DO_V2 -eq 1 ]] && BIN_V2="$(env_bin "$ENV_V2")"

# ---- 3. R packages ----------------------------------------------------------
if [[ $DO_R -eq 1 ]]; then
  step "3/5  R packages (renv::restore())"
  if ! command -v Rscript >/dev/null; then
    warn "R is not installed. Install R 4.6 (https://cran.r-project.org/bin/linux/ubuntu/), then rerun ./install.sh --skip-conda --skip-downloads"
  else
    r_version="$(Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || true)"
    echo "    R $r_version found (renv.lock was recorded with R 4.6)"
    echo "    installing packages; the first time this builds Bioconductor packages and can take 20+ minutes ..."
    if (cd "$REPO_DIR" && Rscript -e 'renv::restore(prompt = FALSE)'); then
      ok "R packages installed"
    else
      warn "renv::restore() failed; see the messages above"
    fi
  fi
fi

# ---- 4. models and reference files -----------------------------------------
REF_DIR="$DATA_DIR/sturgeon_v1_R_scripts"
MODEL_DIR="$DATA_DIR/models"
MODEL_V1="$MODEL_DIR/general.zip"
MODEL_V2=""

if [[ $DO_DOWNLOADS -eq 1 ]]; then
  step "4/5  Sturgeon models and reference files"
  cat <<EOF
    The Sturgeon models and the Sturgeon v1 reference files are not part of
    Shiny-Sturgeon. They are distributed by the Sturgeon authors under their
    own license (an evaluation license for Sturgeon v1), which you can read at
      https://github.com/UMCUGenetics/sturgeon/blob/$STURGEON_V1_COMMIT/LICENSE
      https://github.com/deRidderLab/sturgeon-v2
    By downloading them you agree to those terms.
EOF
  if [[ $ACCEPT_LICENSE -eq 0 ]]; then
    if [[ $ASSUME_YES -eq 1 ]]; then
      fail "with --yes, also pass --accept-sturgeon-license (or use --skip-downloads)"
    fi
    ask "Do you accept the Sturgeon license terms and want to download the files?" \
      || { warn "downloads skipped; fill in the reference files on the Setup tab yourself"; DO_DOWNLOADS=0; }
  fi
fi

if [[ $DO_DOWNLOADS -eq 1 ]]; then
  mkdir -p "$REF_DIR" "$MODEL_DIR"
  download "https://raw.githubusercontent.com/UMCUGenetics/sturgeon/$STURGEON_V1_COMMIT/LICENSE" \
    "$DATA_DIR/STURGEON_LICENSE.txt" || true
  for f in "${STURGEON_V1_FILES[@]}"; do
    download "$STURGEON_V1_RAW/R_scripts/$f" "$REF_DIR/$f" || true
  done

  download "$MODEL_V1_URL" "$MODEL_V1" || true
  if [[ -e "$MODEL_V1" ]] && ! is_zip "$MODEL_V1"; then
    rm -f "$MODEL_V1"
    warn "the Sturgeon v1 model download was not a zip file. Download general.zip by hand from the Sturgeon README and put it at $MODEL_V1"
  fi

  if [[ $DO_V2 -eq 1 ]]; then
    # Ask Zenodo's API for the record's files, so a renamed zip still works
    py="$(command -v python3 || true)"; [[ -z "$py" && -n "$BIN_V2" ]] && py="$BIN_V2/python"
    listing="$(curl -fsSL -A "$USER_AGENT" "https://zenodo.org/api/records/$MODEL_V2_ZENODO_RECORD" 2>/dev/null \
      | "${py:-python3}" -c 'import json,sys
for f in json.load(sys.stdin).get("files", []):
    if f["key"].endswith(".zip"): print(f["key"] + "\t" + f["links"]["self"]); break' 2>/dev/null || true)"
    if [[ -n "$listing" ]]; then
      MODEL_V2="$MODEL_DIR/$(cut -f1 <<<"$listing")"
      download "$(cut -f2 <<<"$listing")" "$MODEL_V2" || MODEL_V2=""
    else
      existing="$(ls "$MODEL_DIR"/cns-v2*.zip 2>/dev/null | head -1 || true)"
      if [[ -n "$existing" ]]; then
        MODEL_V2="$existing"; ok "already present: $MODEL_V2"
      else
        warn "could not reach Zenodo record $MODEL_V2_ZENODO_RECORD (Zenodo sometimes blocks script downloads; a browser usually works). Download the model zip from https://zenodo.org/records/$MODEL_V2_ZENODO_RECORD into $MODEL_DIR and rerun ./install.sh --skip-system --skip-conda --skip-r"
      fi
    fi
  fi
else
  # keep paths to files a previous run already downloaded
  [[ $DO_V2 -eq 1 ]] && MODEL_V2="$(ls "$MODEL_DIR"/cns-v2*.zip 2>/dev/null | head -1 || true)"
fi

# ---- 5. config_local.yml ----------------------------------------------------
step "5/5  Writing config_local.yml"
yq() { local s="${1//\\/\\\\}"; printf '"%s"' "${s//\"/\\\"}"; }   # YAML double-quoted string
existing_file() { [[ -n "$1" && -f "$1" ]] && echo "$1" || echo ""; }
tool_path() {  # tool_path BIN_DIR NAME DEFAULT
  if [[ -n "$1" && -x "$1/$2" ]]; then echo "$1/$2"
  else
    [[ -n "$1" ]] && warn "$2 not found in $1"
    echo "$3"
  fi
}

LOCAL_CFG="$REPO_DIR/config_local.yml"
[[ -f "$LOCAL_CFG" ]] && cp "$LOCAL_CFG" "$LOCAL_CFG.bak"

extra_dirs=()
[[ -n "$BIN_V1" ]] && extra_dirs+=("$BIN_V1")
[[ -n "$BIN_V2" ]] && extra_dirs+=("$BIN_V2")

ref() { existing_file "$REF_DIR/$1"; }
model_v2="$(existing_file "$MODEL_V2")"
{
  echo "# Written by install.sh on $(date '+%Y-%m-%d %H:%M'). Machine-specific; not in git."
  echo "# Loaded on top of config_default.yml. Rerun ./install.sh to regenerate it"
  echo "# (the previous version is kept as config_local.yml.bak)."
  echo "samtools_bin: $(yq "$(tool_path "$BIN_V1" samtools samtools)")"
  echo "modkit_bin: $(yq "$(tool_path "$BIN_V1" modkit modkit)")"
  echo "sturgeon_bin: $(yq "$(tool_path "$BIN_V1" sturgeon sturgeon)")"
  if [[ $DO_V2 -eq 1 ]]; then
    echo "sturgeon_v2_bin: $(yq "$(tool_path "$BIN_V2" sturgeon-v2 sturgeon-v2)")"
  fi
  if [[ ${#extra_dirs[@]} -gt 0 ]]; then
    echo "extra_path_dirs:"
    for d in "${extra_dirs[@]}"; do echo "  - $(yq "$d")"; done
  fi
  echo "probes_file: $(yq "$(ref probelocs_chm13.bed)")"
  echo "model_file_v1: $(yq "$(existing_file "$MODEL_V1")")"
  echo "bins_rds: $(yq "$(ref chm13_v2.0_bins.rds)")"
  echo "cnvplot_reffile: $(yq "$(ref bins_export.bed)")"
  echo "relevant_genes_bed: $(yq "$(ref relevant_genes_with_chm13v2_1Mb_bin_nrs.bed)")"
  echo "color_translation_rds: $(yq "$(ref color_translation.rds)")"
  echo "class_mapping_v2: $(yq "$REPO_DIR/class_mapping_v2.csv")"
  echo "model_file_v2: $(yq "$model_v2")"
} > "$LOCAL_CFG"
ok "$LOCAL_CFG"

for f in probelocs_chm13.bed chm13_v2.0_bins.rds bins_export.bed color_translation.rds; do
  [[ -n "$(ref "$f")" ]] || warn "reference file missing: $REF_DIR/$f"
done
[[ -n "$(existing_file "$MODEL_V1")" ]] || warn "Sturgeon v1 model missing: $MODEL_V1"
[[ $DO_V2 -eq 1 && -z "$model_v2" ]] && warn "no Sturgeon v2 model; the app will run v1 only"

# ---- summary ----------------------------------------------------------------
echo
if [[ ${#WARNINGS[@]} -eq 0 ]]; then
  echo "${GREEN}${BOLD}Installation complete.${RESET}"
else
  echo "${YELLOW}${BOLD}Installation finished with ${#WARNINGS[@]} warning(s):${RESET}"
  for w in "${WARNINGS[@]}"; do echo "  - $w"; done
fi
cat <<EOF

Start the app with:
  $REPO_DIR/run_app.sh

On the Setup tab, fill in the MinKNOW folder, the barcode and an output
folder, then click "Check installation" to confirm everything is found.
EOF
