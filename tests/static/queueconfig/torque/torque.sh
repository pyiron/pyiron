#!/bin/bash
#PBS -q queue@host
#PBS -l ppn={{cores}}
#PBS -N {{job_name}}
{%- if memory_max %}
#PBS -l mem={{memory_max}}
{%- endif %}
{%- if run_time_max %}
#PBS -l walltime={{run_time_max}}
{%- endif %}
#PBS -V
#PBS -d {{working_directory}}

{{command}}