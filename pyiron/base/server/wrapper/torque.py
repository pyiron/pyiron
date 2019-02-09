class TorqueCommands(object):
    @property
    def submit_job_command(self):
        return ['qsub', '-terse']

    @property
    def delete_job_command(self):
        return ['qdel']

    @property
    def enable_reservation_command(self):
        raise NotImplementedError()

    @property
    def get_queue_status_command(self):
        return ['qstat', '-x']

    @staticmethod
    def convert_queue_status(queue_status_output):
        raise NotImplementedError()
