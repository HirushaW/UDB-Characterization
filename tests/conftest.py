import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLE = os.path.join(ROOT, 'examples', 'pamam_g2', 'input')

if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
