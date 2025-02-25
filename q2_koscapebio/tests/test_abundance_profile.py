import subprocess
import os
from qiime2 import Artifact
import pytest
import pandas as pd
import shutil

current_dir = os.path.dirname(__file__)
raw_table_path = os.path.join(current_dir, 'table_final.qza')
qiime_table_path = os.path.join(current_dir, 'table.qza')
expected_result_path = os.path.join(current_dir, 'rel_abundance.qza')

intermediate_output_dir = os.path.join(current_dir, 'rel_abundance_files')
output_dir = os.path.join(current_dir, 'abundance_profile_output')

def test_abundance_profile():
    for dir_path in [intermediate_output_dir, output_dir]:
        if os.path.exists(dir_path):
            shutil.rmtree(dir_path)
        os.makedirs(dir_path, exist_ok=True)

    try:
        cmd = [
            'qiime', 'koscapebio', 'abundance-profile',
            '--p-raw-table', raw_table_path,
            '--p-qiime-table', qiime_table_path,
            '--p-work-dir', intermediate_output_dir,
            '--o-relative-abundance', os.path.join(output_dir, 'test_rel_abundance.qza')
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Command failed: {result.stderr}"
        expected_result = Artifact.load(expected_result_path)
        actual_result = Artifact.load(os.path.join(output_dir, 'test_rel_abundance.qza'))
        expected_df = expected_result.view(pd.DataFrame)
        actual_df = actual_result.view(pd.DataFrame)
        pd.testing.assert_frame_equal(actual_df, expected_df)

    except Exception as e:
        with open(os.path.join(current_dir, 'test_abundance_profile.log'), 'w') as log_file:
            log_file.write(f"Test failed: {e}")
        raise
    else:
        shutil.rmtree(intermediate_output_dir)
        shutil.rmtree(output_dir)

