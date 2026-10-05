import os
import subprocess

root = "results"

for folder in os.listdir(root):
    path = os.path.join(root, folder)

    if not os.path.isdir(path):
        continue

    subprocess.run([
        "sbatch",
        "rotacf_job_template.sbatch",
        path
    ])
