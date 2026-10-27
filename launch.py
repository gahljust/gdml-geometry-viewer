#!/usr/bin/env python3
"""Start the detector workspace."""
import runpy
import sys
from pathlib import Path

if __name__ == '__main__':
    if sys.version_info < (3, 9):
        sys.exit('Python 3.9 or newer is required.')
    sys.dont_write_bytecode = True
    runpy.run_path(str(Path(__file__).resolve().parent / 'geometry/run.py'), run_name='__main__')
