#!/bin/sh
# properties = {"type": "single", "rule": "dataPrepAutomate", "local": false, "input": ["code/dataPrepAutomated.sh"], "output": ["/home1/zhuyixin/zhuyixin_proj/AssmQuality/aligned_bam/test/test_merged_sorted.bam", "/home1/zhuyixin/zhuyixin_proj/AssmQuality/aligned_bam/test/test_merged_sorted.bam.csi"], "wildcards": {"HOME": "/home1/zhuyixin/zhuyixin_proj/AssmQuality", "species": "test"}, "params": {"species": "test", "source": "hifi_fastq", "haploid": "True", "conda": "/spack/conda/miniconda3/23.10.0/etc/profile.d/conda.sh"}, "log": ["/home1/zhuyixin/zhuyixin_proj/AssmQuality/log/dataPrepAutomated.test.log"], "threads": 10, "resources": {"mem_mb": 60000, "mem_mib": 954, "disk_mb": 1000, "disk_mib": 954, "tmpdir": "<TBD>", "mem": "60G"}, "jobid": 3, "cluster": {}}
cd /home1/zhuyixin/ig-assembly-eval && /home1/zhuyixin/.conda/envs/assembly/bin/python3.10 -m snakemake --snakefile '/home1/zhuyixin/ig-assembly-eval/Snakefile' --target-jobs 'dataPrepAutomate:HOME=/home1/zhuyixin/zhuyixin_proj/AssmQuality,species=test' --allowed-rules 'dataPrepAutomate' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=60000' 'mem_mib=954' 'disk_mb=1000' 'disk_mib=954' --wait-for-files '/home1/zhuyixin/ig-assembly-eval/.snakemake/tmp.5xn7134p' 'code/dataPrepAutomated.sh' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'software-env' 'code' 'params' 'input' 'mtime' --skip-script-cleanup  --conda-frontend 'mamba' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --printshellcmds  --latency-wait 60000 --scheduler 'ilp' --scheduler-solver-path '/home1/zhuyixin/.conda/envs/assembly/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && touch '/home1/zhuyixin/ig-assembly-eval/.snakemake/tmp.5xn7134p/3.jobfinished' || (touch '/home1/zhuyixin/ig-assembly-eval/.snakemake/tmp.5xn7134p/3.jobfailed'; exit 1)

