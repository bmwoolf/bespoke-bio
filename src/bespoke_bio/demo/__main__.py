import argparse
import json
from pathlib import Path

from .catalog import build_index

parser = argparse.ArgumentParser(description="Index a verified Stage 1 package for the demo")
parser.add_argument("package", type=Path)
parser.add_argument("index", type=Path)
arguments = parser.parse_args()
print(json.dumps(build_index(arguments.package, arguments.index), indent=2))
