#!/bin/bash
#SBATCH -p batch                                                # partition (this is the queue your job will be added to) 
#SBATCH -N 1                                                    # number of nodes (use a single node)
#SBATCH -n 8                                                    # number of cores (sequential job => uses 1 core)
#SBATCH --time=24:00:00                                         # time allocation, which has the format (D-HH:MM:SS), here set to 1 hour
#SBATCH --mem=16GB 

bash GetMapReads.sh
bash CreateTSV.sh
python TSVlength.py
python main.py 
