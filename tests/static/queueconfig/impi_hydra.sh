#!/bin/bash
#$ -N {{job_name}}
#$ -wd {{working_directory}}
{%- if cores %}
#$ -pe impi_hydra {{cores}}
{%- endif %}
{%- if memory_max %}
#$ -l h_vmem={{memory_max}}
{%- endif %}
{%- if run_time_max %}
#$ -l h_rt={{run_time_max}}
{%- endif %}
#$ -o time.out
#$ -e error.out

{{command}}
