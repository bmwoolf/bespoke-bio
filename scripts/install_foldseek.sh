#!/usr/bin/env bash
set -euo pipefail
mkdir -p vendor
case "$(uname -s)" in
  Darwin) url=https://mmseqs.com/foldseek/foldseek-osx-universal.tar.gz ;;
  Linux) url=https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz ;;
  *) echo 'Install Foldseek manually on this platform'; exit 1 ;;
esac
curl --fail --location "$url" --output vendor/foldseek.tar.gz
tar -xzf vendor/foldseek.tar.gz -C vendor
vendor/foldseek/bin/foldseek version
# This URL is mutable. Record the printed version and archive hash in every deployment.
