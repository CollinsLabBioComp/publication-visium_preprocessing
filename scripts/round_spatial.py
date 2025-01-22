import numpy as np
import pandas as pd
import math
import argparse
import os


parser = argparse.ArgumentParser(description='Round spatial')
parser.add_argument('--dataDir', type=str, nargs=1, default='/home/'+os.environ['USER']+'/DeepSpaCE/data', required=True,
                    help='Data directory (default: '+'/home/'+os.environ['USER']+'/DeepSpaCE/data'+')')

args = parser.parse_args()
dataDir = args.dataDir[0]


tissue = pd.read_csv(dataDir+'spatial/tissue_positions.csv', header=None)
tissue[4] = tissue[4].round().astype(int)
tissue[5] = tissue[5].round().astype(int)
tissue.to_csv(dataDir+'spatial/tissue_positions.csv', index=False, header=None)

