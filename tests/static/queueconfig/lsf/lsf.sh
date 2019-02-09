#!/bin/bash
#BSUB -q queue
#BSUB -J {{job_name}}
#BSUB -o time.out
#BSUB -n {{cores}}
#BSUB -cwd {{working_directory}}
#BSUB -e error.out
{%- if run_time_max %}
#BSUB -W {{run_time_max}}
{%- endif %}
{%- if memory_max %}
#BSUB -M {{memory_max}}
{%- endif %}

{{command}}