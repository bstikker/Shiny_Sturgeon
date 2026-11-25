#!/usr/bin/env bash
# Start Shiny-Sturgeon and open it in the web browser.
#   ./run_app.sh            # on http://127.0.0.1:8080
#   ./run_app.sh 9000       # on another port
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
PORT="${1:-8080}"
command -v Rscript >/dev/null || { echo "R is not installed; see docs/install-guide.md"; exit 1; }
echo "Starting Shiny-Sturgeon on http://127.0.0.1:$PORT  (press Ctrl+C here to stop it)"
exec Rscript -e "shiny::runApp('.', host = '127.0.0.1', port = $PORT, launch.browser = interactive() || nzchar(Sys.getenv('DISPLAY')))"
