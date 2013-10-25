import sys, os
cur_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(cur_dir, 'python'))
import pyHyp
__all__ = ['pyHyp']
