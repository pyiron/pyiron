# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.
import os
import h5py
import shutil
from pyiron import Project
from pyiron_base.settings.generic import Settings


s = Settings()

def register(parser):
    parser.add_argument(
            "-i", "--input-path", type = os.path.abspath,
    )
    parser.add_argument(
            "-o", "--output-path", type = os.path.abspath,
    )

def main(args):
    with h5py.File(args.input_path, mode="r") as f:
        job_name = list(f.keys())[0]
    project_path = os.path.join(os.path.abspath("."), job_name + ".h5")
    shutil.copy(args.input_path, project_path)

    file = os.path.basename(project_path)
    job_name = os.path.splitext(file)[0]

    db_project_path = s.top_path(project_path)
    project = os.path.dirname(project_path)
    db_project = (project + "/").replace(db_project_path, "")
    job_reload = Project(project).load_from_jobpath(
        job_id=None,
        db_entry={
            "id": 1000,
            "status": "",
            "chemicalformula": "",
            "job": job_name,
            "subjob": "/" + job_name,
            "projectpath": db_project_path,
            "project": db_project,
            "hamilton": "",
            "hamversion": "",
            "parentid": None,
            "masterid": None,
        },
        convert_to_object=True,
    )
    job_reload.status.initialized = True
    job_reload.server.run_mode.modal = True
    job_reload.run()
    shutil.copy(project_path, args.output_path)
