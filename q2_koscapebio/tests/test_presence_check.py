import subprocess
import os
from qiime2 import Artifact
import pytest
import pandas as pd
import shutil

current_dir = os.path.dirname(__file__)
rep_seq_path = os.path.join(current_dir, 'rep-seqs.qza')
table_path = os.path.join(current_dir, 'table.qza')
expected_result_path = os.path.join(current_dir, 'table_final.qza')

output_dir = os.path.join(current_dir, 'presence_check_output')

def test_presence_check():
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    try:
        cmd = [
            'qiime', 'koscapebio', 'presence-check',
            '--p-rep-seq-path', rep_seq_path,
            '--p-table-path', table_path,
            '--p-perc-identity', '1',
            '--p-strand','both',
            '--output-dir', output_dir
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Command failed: {result.stderr}"
        
        expected_result = Artifact.load(expected_result_path)
        actual_result = Artifact.load(os.path.join(output_dir, 'clustered_table.qza'))
        
        expected_df = expected_result.view(pd.DataFrame)
        actual_df = actual_result.view(pd.DataFrame)
        
        pd.testing.assert_frame_equal(actual_df, expected_df)

    except Exception as e:
        with open(os.path.join(current_dir, 'test_presence_check.log'), 'w') as log_file:
            log_file.write(f"Test failed: {e}")
        raise
    else:
        shutil.rmtree(output_dir)

