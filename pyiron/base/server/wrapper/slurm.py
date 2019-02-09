class SlurmCommands(object):
    @property
    def submit_job_command(self):
        return ['sbatch', '--parsable']

    @property
    def delete_job_command(self):
        return ['scancel']

    @property
    def enable_reservation_command(self):
        raise NotImplementedError()

    @property
    def get_queue_status_command(self):
        return ['squeue']

    @staticmethod
    def convert_queue_status(queue_status_output):
        raise NotImplementedError()
