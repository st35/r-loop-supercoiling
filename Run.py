import sys
from os.path import isfile
import subprocess
from pathlib import Path

foldername_0 = 'outputfiles'
foldername_1 = 'inputfiles'

if not Path(foldername_0).is_dir():
	Path(foldername_0).mkdir()

if not Path(foldername_1).is_dir():
	Path(foldername_1).mkdir()

restart_flag = 0

start_rank = 16

import numpy as np

promoter_flag = [0, 1, 2]
GQ_on_val = [1.0, 6.0, 10.0]
eR_on_val = [0.1, 0.6, 1.0, 6.0]
topo_val = [10.0, 60.0, 100.0]
PQS_flag = [0, 1]

for plasmid_flag in range(2):
	for promoter_id in range(3):
		for GQ_on_id in range(len(GQ_on_val)):
			for eR_on_id in range(len(eR_on_val)):
				for topo_id in range(len(topo_val)):
					for PQS_id in range(2):
						folder_0 = 'outputfiles/RUN_' + str(plasmid_flag) + '_' + str(promoter_flag[promoter_id]) + '_' + str(GQ_on_id) + '_' + str(eR_on_id) + '_' + str(topo_id) + '_' + str(PQS_id)
						folder_1 = 'inputfiles/RUN_' + str(plasmid_flag) + '_' + str(promoter_flag[promoter_id]) + '_' + str(GQ_on_id) + '_' + str(eR_on_id) + '_' + str(topo_id) + '_' + str(PQS_id)
						if not Path(folder_0).is_dir():
							Path(folder_0).mkdir()
						if not Path(folder_1).is_dir():
							Path(folder_1).mkdir()
						marker_filename = 'outputfiles/RUN_' + str(plasmid_flag) + '_' + str(promoter_flag) + '_' + str(GQ_on_id) + '_' + str(eR_on_id) + '_' + str(topo_id) + '_' + str(PQS_flag) + '/rates_0.log'
						if isfile(marker_filename):
							continue
						A = [str(restart_flag), str(plasmid_flag), str(promoter_flag[promoter_id]), str(GQ_on_id), str(GQ_on_val[GQ_on_id]), str(eR_on_id), str(eR_on_val[eR_on_id]), str(topo_id), str(topo_val[topo_id]), str(PQS_flag[PQS_id])]
						print(A)
						slurm_out = subprocess.run(['sbatch', 'submit_job.slurm', 'configs/test.config', str(restart_flag), str(plasmid_flag), str(promoter_flag[promoter_id]), str(GQ_on_id), str(GQ_on_val[GQ_on_id]), str(eR_on_id), str(eR_on_val[eR_on_id]), str(topo_id), str(topo_val[topo_id]), str(PQS_flag[PQS_id]), str(start_rank)], capture_output = True, text = True)
						print(slurm_out.stdout, end = '')
