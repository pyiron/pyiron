import os
import pandas
import unittest
from pyiron.base.server.queueadapter import QueueAdapter


class TestRunmode(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path = os.path.dirname(os.path.abspath(__file__))
        cls.torque = QueueAdapter(directory=os.path.join(cls.path, '../../static/queueconfig/torque'))
        cls.slurm = QueueAdapter(directory=os.path.join(cls.path, '../../static/queueconfig/slurm'))
        cls.lsf = QueueAdapter(directory=os.path.join(cls.path, '../../static/queueconfig/lsf'))
        cls.sge = QueueAdapter(directory=os.path.join(cls.path, '../../static/queueconfig/sge'))
        cls.moab = QueueAdapter(directory=os.path.join(cls.path, '../../static/queueconfig/moab'))

    def test_config(self):
        self.assertEqual(self.torque.config['queue_type'], 'TORQUE')
        self.assertEqual(self.slurm.config['queue_type'], 'SLURM')
        self.assertEqual(self.lsf.config['queue_type'], 'LSF')
        self.assertEqual(self.sge.config['queue_type'], 'SGE')
        self.assertEqual(self.moab.config['queue_type'], 'MOAB')
        self.assertEqual(self.torque.config['queue_primary'], 'torque')
        self.assertEqual(self.slurm.config['queue_primary'], 'slurm')
        self.assertEqual(self.lsf.config['queue_primary'], 'lsf')
        self.assertEqual(self.sge.config['queue_primary'], 'impi_hydra_small')
        self.assertEqual(self.moab.config['queue_primary'], 'moab')

    def test_value_in_range(self):
        self.assertEqual(None, self.sge._value_in_range(value=None, value_min=None, value_max=None))
        self.assertEqual(1, self.sge._value_in_range(value=None, value_min=1, value_max=None))
        self.assertEqual(1, self.sge._value_in_range(value=None, value_min=None, value_max=1))
        self.assertEqual(1, self.sge._value_in_range(value=1, value_min=None, value_max=None))
        self.assertEqual(1, self.sge._value_in_range(value=0, value_min=1, value_max=None))
        self.assertEqual(1, self.sge._value_in_range(value=2, value_min=None, value_max=1))
        self.assertEqual(1, self.sge._value_in_range(value=1, value_min=0, value_max=2))

    def test_job_submission_template(self):
        self.assertRaises(ValueError, self.sge._job_submission_template, command=None)
        self.assertRaises(TypeError, self.sge._job_submission_template, command=1)
        template = "#!/bin/bash\n#$ -N job.py\n#$ -wd .\n#$ -pe impi_hydra_small 1\n#$ -l h_rt=604800\n" \
                   "#$ -o time.out\n#$ -e error.out\n\npython test.py"
        self.assertEqual(self.sge._job_submission_template(command='python test.py'), template)
        template = "#!/bin/bash\n#BSUB -q queue\n#BSUB -J job.py\n#BSUB -o time.out\n#BSUB -n 10\n#BSUB -cwd .\n" \
                   "#BSUB -e error.out\n#BSUB -W 259200\n\npython test.py"
        self.assertEqual(self.lsf._job_submission_template(command='python test.py'), template)

    def test_interfaces(self):
        self.assertEqual(self.sge._commands.submit_job_command, ['qsub', '-terse'])
        self.assertEqual(self.sge._commands.delete_job_command, ['qdel'])
        self.assertEqual(self.sge._commands.enable_reservation_command, ['qalter', '-R', 'y'])
        self.assertEqual(self.sge._commands.get_queue_status_command, ['qstat', '-xml'])
        self.assertEqual(self.torque._commands.submit_job_command, ['qsub', '-terse'])
        self.assertEqual(self.torque._commands.delete_job_command, ['qdel'])
        self.assertEqual(self.torque._commands.get_queue_status_command, ['qstat', '-x'])
        self.assertEqual(self.lsf._commands.submit_job_command, ['bsub', '-terse'])
        self.assertEqual(self.lsf._commands.delete_job_command, ['bkill'])
        self.assertEqual(self.lsf._commands.get_queue_status_command, ['qstat', '-x'])
        self.assertEqual(self.slurm._commands.submit_job_command, ['sbatch', '--parsable'])
        self.assertEqual(self.slurm._commands.delete_job_command, ['scancel'])
        self.assertEqual(self.slurm._commands.get_queue_status_command, ['squeue'])
        self.assertEqual(self.moab._commands.submit_job_command, ['msub'])
        self.assertEqual(self.moab._commands.delete_job_command, ['mjobctl', '-c'])
        self.assertEqual(self.moab._commands.get_queue_status_command, ['mdiag', '-x'])

    def test_convert_queue_status(self):
        with open(os.path.join(self.path, '../../static/queueconfig/sge', 'qstat.xml'), 'r') as f:
            content = f.read()
        df_running = pandas.DataFrame({'jobid': ['2836045'], 'user': ['friko'], 'jobname': ['vasp.5.3.5'],
                                       'status': ['running']})
        df_pending = pandas.DataFrame({'jobid': ['2836046', '2967274'], 'user': ['friko', 'janj'],
                                       'jobname': ['vasp.5.3.5', 'hello.py'], 'status': ['pending', 'error']})
        df_merge = df_running.append(df_pending, sort=True)
        df = pandas.DataFrame({'jobid': pandas.to_numeric(df_merge.jobid), 'user': df_merge.user,
                               'jobname': df_merge.jobname, 'status': df_merge.status})
        self.assertTrue(df.equals(self.sge._commands.convert_queue_status(queue_status_output=content)))
